# Read in adult_aged_bodymap data
load('adult_aged_bodymap/bodymap_links.RData')
# Read in TAC data
load('F1_TAC_Sarah/TAC_SHAM_links.RData')
# Overlap lists
load('LINKS_study/adult_aged_TAC_link_overlaps.RData')



tmp <- merge(bodymap[['He_78w']], TAC_SHAM[['He_TAC28d_RNA_XistBxC_-+_XX']], by.x = c('name_base', 'name_target'), by.y = c('name_base', 'name_target'))
delta_AR <- abs(tmp$allelic_ratio_base.x - tmp$allelic_ratio_base.y)
names(delta_AR) <- tmp$name_base
lapply(overlap_lists[['aged_tac_only_enh']], function(x){
    base <- strsplit(x, '\\|')[[1]][1]
    delta_AR[base]
})


9W atrialCM
9w FB
78w SMC
78w atrialCM
78w LEC
78w EC