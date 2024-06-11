def classify_sample(sample, files) {

    def meta = [:]
    meta.is_paired = (files instanceof Collection && files.size() == 2)
    meta.id = sample

    return [meta, files]

    if (meta.is_paired) {
        return [meta, files]
    }

    return [meta, [files]]

}


def classify_sample_with_library_info(sample, library, files) {
    
    def meta = [:]
    meta.is_paired = (files instanceof Collection && files.size() == 2)
    meta.id = sample
    meta.library = (library == "1") ? "paired" : "single"

    return [meta, files]

}