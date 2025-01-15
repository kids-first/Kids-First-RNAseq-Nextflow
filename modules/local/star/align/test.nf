include { STAR_ALIGN } from './main.nf'

workflow {
    main:
    genomeDir = Channel.fromPath(params.genomeDir)
    readFilesCommand = Channel.value(params.readFilesCommand)
    readFilesManifest = Channel.fromPath(params.readFilesManifest)
    outFileNamePrefix = Channel.value(params.outFileNamePrefix)
    STAR_ALIGN(genomeDir, readFilesCommand, readFilesManifest, outFileNamePrefix)
}