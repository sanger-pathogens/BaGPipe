def validate_parameters() {
    def errors = 0

    if (params.manifest) {
        manifest=file(params.manifest)
        if (!manifest.exists()) {
            log.error("The genome manifest file specified does not exist.")
            errors += 1
        }
    } else {
        log.error("No genome manifest file specified. Please specify one using the --manifest option.")
        errors += 1
    }

    if (!["unitig", "pa", "snp"].contains(params.genotype_method)) {
        log.error("Invalid genotype method. Please use one of the three options: unitig|pa|snp .")
        errors += 1
    }

    if (params.reference) {
        reference=file(params.reference)
        if (!reference.exists()) {
            log.error("The reference manifest file specified does not exist.")
            errors += 1
        }
    }

    if (params.mygff) {
        gff=file(params.mygff)
        if (!gff.exists()) {
            log.error("The annotated genome (gff) manifest file specified does not exist.")
            errors += 1
        }
    }

    if (errors > 0) {
        log.error(String.format("%d errors detected", errors))
        exit 1
    }
}
