version 1.0

workflow callRemora {

    input {
        Array[File] inputTarballOrFast5s
        File referenceFasta

        # Guppy configuration
        String guppyConfig = "dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg"

        # naming for final output
        String? sampleIdentifier

        # resources for guppy
        Int memSizeGB = 128
        Int threadCount = 64
        Int gpuCount = 0
        Array[String] zones = ['us-central1-c']
        String dockerImage = "xiaoyuz/guppy:latest"
    }

    scatter (inputFile in inputTarballOrFast5s) {
        call untar {
            input:
                fileToUntar=inputFile,
                zones=zones,
                dockerImage=dockerImage
        }

        Array[File] fast5s = if untar.didUntar then untar.untarredFast5s else [inputFile]

        if (gpuCount > 0) {
            call remoraGPU {
                input:
                    inputFast5s = fast5s,
                    referenceFasta = referenceFasta,
                    guppyConfig = guppyConfig,
                    memSizeGB = memSizeGB,
                    threadCount = threadCount,
                    diskSizeGB = untar.fileSizeGB * 2 + 5,
                    gpuCount = gpuCount,
                    zones=zones,
                    dockerImage=dockerImage
            }
        }

        if (gpuCount <= 0) {
            call remoraCPU {
                input:
                    inputFast5s = fast5s,
                    referenceFasta = referenceFasta,
                    guppyConfig = guppyConfig,
                    memSizeGB = memSizeGB,
                    threadCount = threadCount,
                    diskSizeGB = untar.fileSizeGB * 2 + 5,
                    zones=zones,
                    dockerImage=dockerImage
            }
        }

    }

    call sum {
        input:
            integers = if (gpuCount > 0) then remoraGPU.fileSizeGB else remoraCPU.fileSizeGB,
            dockerImage=dockerImage
    }

    call mergeRemora {
        input:
            sampleIdentifier = sampleIdentifier,
            remoraOutputTarballs = if (gpuCount > 0) then remoraGPU.outputTarball else remoraCPU.outputTarball,
            diskSizeGB = sum.value * 5, #output tar, output untar, merged files, tarred merge, slop
            zones=zones,
            dockerImage=dockerImage
    }

    output {
        File mergedRemoraResults = mergeRemora.mergedTarball
    }
}

task untar {
    input {
        File fileToUntar
        Int diskSizeGB = 512
        Array[String] zones = ['us-central1-c']
        String dockerImage = "xiaoyuz/guppy:latest"
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # untar or don't
        mkdir tmp
        cd tmp
        if [[ "~{fileToUntar}" == *.tar ]] || [[ "~{fileToUntar}" == *.tar.gz ]] ; then
            tar xvf ~{fileToUntar}
            echo "true" >../untarred
        else
            echo "false" >../untarred
        fi
        cd ..

        # move everything to output
        mkdir output
        for FILE in `find tmp/ -name "*.fast5"` ; do
            mv $FILE output
        done

        # get output size
        if [[ `ls output | wc -l` == 0 ]] ; then
            OUTPUTSIZE=`du -s -BG ~{fileToUntar} | sed 's/G.*//'`
        else
            OUTPUTSIZE=`du -s -BG output/ | sed 's/G.*//'`
        fi
        echo $OUTPUTSIZE >outputsize
    >>>

    output {
        Boolean didUntar = read_boolean("untarred")
        Array[File] untarredFast5s = glob("output/*")
        Int fileSizeGB = read_int("outputsize")
    }

    runtime {
        memory: "2 GB"
        cpu: 1
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 1
        zones: zones
    }

}

task remoraGPU {
    input {
        # files
        Array[File] inputFast5s
        File referenceFasta

        # remora configuration

        String guppyConfig = "dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg"

        # resources
        Int memSizeGB = 64
        Int threadCount = 12
        Int diskSizeGB = 128
        Int gpuCount = 1
        String gpuType = "nvidia-tesla-v100"
        String nvidiaDriverVersion = "418.87.00"
        Int maxRetries = 4 # workaround for Terra failure to initilize drivers
        Array[String] zones = 	[ "us-central1-c" ]
        String dockerImage = "xiaoyuz/guppy:latest"
    }

	command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # prepare input
        mkdir input_files
        for FILE in ~{ sep=' ' inputFast5s } ; do
            ln -s $FILE input_files/
        done

        # logging
        ls -lahR
        df -H

        # start constructing remora command
        cmd=(/ont-guppy/bin/guppy_basecaller -i input_files/)
        cmd+=( -a ~{referenceFasta} )
        cmd+=( -s output/ )
        cmd+=( -c dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg )
        cmd+=( --bam_out )

        # add GPU numbers
        cmd+=( --device cuda:all)

        # run remora command
        "${cmd[@]}"

        # save output
        UUID=`uuid`
        mkdir output_$UUID
        ls output/ | xargs -n1 -I{} mv output/{} output_$UUID/{}
        tar czvf remora_output_$UUID.tar.gz output_$UUID/

        # get output size
        du -s -BG output_$UUID/ | sed 's/G.*//' >outputsize

	>>>
	output {
		File outputTarball = glob("remora_output_*.tar.gz")[0]
        Int fileSizeGB = read_int("outputsize")
	}
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        gpuCount: gpuCount
        gpuType: gpuType
        maxRetries: maxRetries
        nvidiaDriverVersion: nvidiaDriverVersion
        docker: dockerImage
        zones: zones
    }
}

task remoraCPU {
    input {
        # files
        Array[File] inputFast5s
        File referenceFasta

        # remora configuration

        String guppyConfig = "dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg"

        # resources
        Int memSizeGB = 64
        Int threadCount = 12
        Int diskSizeGB = 128
        Int gpuCount = 1
        Int maxRetries = 4 # workaround for Terra failure to initilize drivers
        Array[String] zones =   [ "us-central1-c" ]
        String dockerImage = "xiaoyuz/guppy:latest"
    }

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # prepare input
        mkdir input_files
        for FILE in ~{ sep=' ' inputFast5s } ; do
            ln -s $FILE input_files/
        done

        # logging
        ls -lahR
        df -H

        # start constructing remora command
        cmd=(/ont-guppy/bin/guppy_basecaller -i input_files/)
        cmd+=( -a ~{referenceFasta} )
        cmd+=( -s output/ )
        cmd+=( -c dna_r9.4.1_450bps_modbases_5mc_cg_sup_prom.cfg )
        cmd+=( --bam_out )

        # run remora command
        echo "$cmd"
        "${cmd[@]}"

        # save output
        UUID=`uuid`
        mkdir output_$UUID
        ls output/ | xargs -n1 -I{} mv output/{} output_$UUID/{}
        tar czvf remora_output_$UUID.tar.gz output_$UUID/

        # get output size
        du -s -BG output_$UUID/ | sed 's/G.*//' >outputsize

    >>>
    output {
        File outputTarball = glob("remora_output_*.tar.gz")[0]
        Int fileSizeGB = read_int("outputsize")
    }
    runtime {
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        maxRetries: maxRetries
        docker: dockerImage
        zones: zones
    }
}

task sum {
    input {
        Array[Int?] integers
        String dockerImage
    }
    Array[Int] allValidIntegers = select_all(integers)

    command <<<
        echo $((0 + ~{sep="+" allValidIntegers}))
    >>>

    output {
        Int value = read_int(stdout())
    }

    runtime {
        preemptible: 1
        docker: dockerImage
    }
}


task mergeRemora {
    input {
        Array[File?] remoraOutputTarballs
        String? sampleIdentifier
        Int threadCount = 8
        Int memSizeGB = 8
        Int diskSizeGB = 128
        Array[String] zones = ['us-central1-c']
        String dockerImage = "xiaoyuz/guppy:latest"
    }
    Array[File] allValidRemoraOutputTarballs = select_all(remoraOutputTarballs)

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # extract tarballs
        mkdir extracted
        cd extracted
        for FILE in ~{sep=' ' allValidRemoraOutputTarballs} ; do
            tar xvf $FILE
        done
        cd ..

        # setup output
        mkdir output

        # handle each output type from remora

        echo "BASECALLS: $(date)"
        mkdir -p tmp_basecalls/pass
        mkdir -p tmp_basecalls/fail
        find extracted/*/pass -name *.fastq | xargs -n1 -I{} mv {} tmp_basecalls/pass
        cat tmp_basecalls/pass/*.fastq >output/pass.merged_basecalls.fastq
        find extracted/*/fail -name *.fastq | xargs -n1 -I{} mv {} tmp_basecalls/fail
        cat tmp_basecalls/fail/*.fastq >output/fail.merged_basecalls.fastq

        echo "MAPPINGS: $(date)"
        mkdir -p tmp_mappings/pass
        mkdir -p tmp_mappings/fail
        find extracted/*/pass -name *.bam | xargs -n1 -I{} bash -c 'samtools sort -@~{threadCount} {} >tmp_mappings/pass/$(basename {})'
        samtools merge -@~{threadCount} output/pass.merged_mappings.bam tmp_mappings/pass/*
        samtools index -@~{threadCount} output/pass.merged_mappings.bam
        find extracted/*/fail -name *.bam | xargs -n1 -I{} bash -c 'samtools sort -@~{threadCount} {} >tmp_mappings/fail/$(basename {})'
        samtools merge -@~{threadCount} output/fail.merged_mappings.bam tmp_mappings/fail/*
        samtools index -@~{threadCount} output/fail.merged_mappings.bam

        mkdir tmp_mapping_summary
        #TODO merge mapping summary

        # finalize output
        ID="~{if defined(sampleIdentifier) then sampleIdentifier + "." else ""}"
        cd output/
        ls * | xargs -n1 -I{} mv {} ${ID}{}
        tar czvf ${ID}merged_remora.tar.gz *
        mv *merged_remora.tar.gz ..
    >>>

    output {
        File mergedTarball = glob("*merged_remora.tar.gz")[0]
    }


    runtime {
        memory: memSizeGB + "GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        docker: dockerImage
        preemptible: 0
        zones: zones
    }


}