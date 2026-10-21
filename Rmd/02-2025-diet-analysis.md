02-2025-diet-analysis
================
diana baetscher
2026-10-15

Following the initial bioinformatic analysis in
`01-2025-mifish-analysis.Rmd`

This is a clean Rmd that doesn’t require moving between R and terminal.

``` r
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(stringr)
library(ggplot2)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ readr     2.1.4
    ## ✔ lubridate 1.9.5     ✔ tibble    3.3.1
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
# for adding distribution data
library(rfishbase)
```

``` r
taxonomy <-read_tsv("../Rmd/blast_outputs/NFS_diet_mifish_shorttaxlineage.txt", col_names = F)
```

    ## Rows: 42102 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (6): X1, X2, X13, X14, X16, X17
    ## dbl (11): X3, X4, X5, X6, X7, X8, X9, X10, X11, X12, X15
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# first clean up the header
taxon_df <- taxonomy %>%
  select(-X2, -(X4:X16)) %>%  #remove unnecessary columns
  group_by(X1, X17) %>% # group by the sequence key and the full taxonomy to reduce duplicate entries
  unique() %>% # doing that reduced the number of entries significantly
  rename(ASV=X1, perc_id=X3, taxonomy=X17) %>% #rename headers
  filter(!str_detect(taxonomy, "environmental")) %>% # filter out any environmental samples
  filter(!str_detect(taxonomy, "synthetic")) %>% # filter out any synthetic "samples"
  filter(!str_detect(taxonomy, "bacterium")) %>%
  filter(!str_detect(taxonomy, "bacteria")) %>%
  filter(!str_detect(taxonomy, "cultured")) %>% # formatting the taxonomy variables
  filter(str_detect(taxonomy, ";")) %>%
  separate(taxonomy, into=c("kingdom", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>% 
  filter(!str_detect(species, " x ")) %>% # remove hybrids
   # remove other stuff that isn't on-target or out of range
   filter(!species %in% c("Homo sapiens", "Canis lupus", "Pan paniscus")) %>% # remove species outside of their range
  # filter(!str_detect(species, " sp. ")) %>%
   filter(!genus %in% c("Sus", "Bos", "Gorilla")) #%>% # remove non-marine mammals
  # filter(!order %in% c("Primates"))
```

## what non-fish are here?

``` r
not_Actinopteri <- taxon_df %>%
  ungroup() %>%
  filter(class != "Actinopteri") %>%
  select(species, genus, family, order, class, phylum) %>%
  unique() %>%
  select(species)

not_Actinopteri %>%
  left_join(., taxon_df)
```

    ## Joining with `by = join_by(species)`

    ## # A tibble: 56 × 9
    ##    species             ASV   perc_id kingdom phylum   class   order family genus
    ##    <chr>               <chr>   <dbl> <chr>   <chr>    <chr>   <chr> <chr>  <chr>
    ##  1 Callorhinus ursinus ASV1    100   ""      Chordata Mammal… Carn… Otari… Call…
    ##  2 Callorhinus ursinus ASV1     99.4 ""      Chordata Mammal… Carn… Otari… Call…
    ##  3 Callorhinus ursinus ASV4     98.8 ""      Chordata Mammal… Carn… Otari… Call…
    ##  4 Callorhinus ursinus ASV4     98.2 ""      Chordata Mammal… Carn… Otari… Call…
    ##  5 Callorhinus ursinus ASV5     99.4 ""      Chordata Mammal… Carn… Otari… Call…
    ##  6 Callorhinus ursinus ASV5     98.8 ""      Chordata Mammal… Carn… Otari… Call…
    ##  7 Callorhinus ursinus ASV6    100   ""      Chordata Mammal… Carn… Otari… Call…
    ##  8 Callorhinus ursinus ASV6     99.4 ""      Chordata Mammal… Carn… Otari… Call…
    ##  9 Callorhinus ursinus ASV9     99.4 ""      Chordata Mammal… Carn… Otari… Call…
    ## 10 Callorhinus ursinus ASV9     98.8 ""      Chordata Mammal… Carn… Otari… Call…
    ## # ℹ 46 more rows

## now look at the fish and figure out what taxa are in/out of our range

``` r
to_check_range <- taxon_df %>%
  ungroup() %>%
  filter(class == "Actinopteri") %>%
  select(species, genus, family, order, class, phylum) %>%
  unique()
```

## check ranges for species using rfishbase

``` r
#first, validate species names using rfishbase synonyms
spp_df <- synonyms(to_check_range$species) # this takes a few minutes
```

    ## Joining with `by = join_by(Subfamily, GenCode, FamCode)`
    ## Joining with `by = join_by(FamCode)`
    ## Joining with `by = join_by(Order, Ordnum, Class, ClassNum)`
    ## Joining with `by = join_by(Class, ClassNum)`

``` r
syn <- spp_df %>%
  filter(Status == "synonym")

# this performs the by-row substitution for the synonyms
to_check_range2 <- to_check_range %>%
  left_join(., spp_df, by = c("species" = "synonym")) %>%
  #filter(Status == "synonym") %>%
  mutate(validated_species = ifelse(Status == "synonym", Species, species))

to_check_range2 %>%
  filter(validated_species != species) %>%
  select(species, Species, validated_species)
```

    ## # A tibble: 6 × 3
    ##   species                  Species                 validated_species      
    ##   <chr>                    <chr>                   <chr>                  
    ## 1 Gadus ogac               Gadus macrocephalus     Gadus macrocephalus    
    ## 2 Sebastes matsubarae      Sebastes matsubarai     Sebastes matsubarai    
    ## 3 Cleisthenes herzensteini Cleisthenes pinetorum   Cleisthenes pinetorum  
    ## 4 Kareius bicoloratus      Platichthys bicoloratus Platichthys bicoloratus
    ## 5 Embassichthys bathybius  Microstomus bathybius   Microstomus bathybius  
    ## 6 Lycodopsis pacificus     Lycodes pacificus       Lycodes pacificus

``` r
# this loads a large table into memory
# but doesn't take long
fishbase_spp <- fb_tbl("species")
head(fishbase_spp)
```

    ## # A tibble: 6 × 102
    ##   SpecCode Genus             Species SpeciesRefNo Author FBname PicPreferredName
    ##      <int> <chr>             <chr>          <int> <chr>  <chr>  <chr>           
    ## 1    64588 Aapticheilichthys webste…        78622 (Hube… <NA>   Prweb_u0.jpg    
    ## 2    16239 Aaptosyax         grypus         10431 Rainb… Giant… Aagry_u0.gif    
    ## 3     2347 Abactochromis     labros…        85491 (Trew… <NA>   Melab_u0.jpg    
    ## 4    62612 Abalistes         filame…        54835 Matsu… <NA>   Abfil_u0.jpg    
    ## 5        9 Abalistes         stella…         9770 (Anon… Starr… Abste_u4.jpg    
    ## 6    11144 Abantennarius     analis          6773 Schul… Tailj… Anana_u2.jpg    
    ## # ℹ 95 more variables: PicPreferredNameM <chr>, PicPreferredNameF <chr>,
    ## #   PicPreferredNameJ <chr>, FamCode <int>, Subfamily <chr>, GenCode <int>,
    ## #   SubGenCode <int>, BodyShapeI <chr>, Source <chr>, AuthorRef <int>,
    ## #   Remark <chr>, TaxIssue <int>, Fresh <int>, Brack <int>, Saltwater <int>,
    ## #   DemersPelag <chr>, AirBreathing <chr>, AirBreathingRef <int>, AnaCat <chr>,
    ## #   MigratRef <int>, DepthRangeShallow <int>, DepthRangeDeep <int>,
    ## #   DepthRangeRef <int>, DepthRangeComShallow <int>, DepthRangeComDeep <int>, …

``` r
# get distribution info
spp_distribution <- faoareas(to_check_range2$validated_species) %>%
  select(Species, FAO) %>%
  unique()
```

    ## Joining with `by = join_by(SpecCode)`

``` r
# add column to designate if we will consider a species as "in range"- for this study, this will be NE Pacific and Arctic Ocean
spp_distribution_range <- spp_distribution %>%
  mutate(in_range = ifelse(is.na(FAO), NA, "no"),
         in_range = ifelse(FAO == "Pacific, Northeast", "yes", in_range))
          #in_range = ifelse(FAO == "Arctic Ocean", "yes", in_range)) # not expecting arctic species, but leaving this in

# confirm that the species look appropriate
# add back synonyms that could be in the blastn output
syn_to_addback <- spp_distribution_range %>%
  inner_join(., syn) %>%
  filter(in_range == "yes") %>%
  select(Species, synonym, FAO, in_range) %>%
  unique() %>%
  rename(species = Species)
```

    ## Joining with `by = join_by(Species)`

``` r
# keep just a list of spp names and yes/no/NA for "in range"  - this way we can keep track of what spp didn't have any reference information in fishbase to determine range
spp_in_range <- spp_distribution_range %>%
  select(Species, FAO, in_range) %>%
  filter(in_range == "yes") %>%
  unique() %>%
  rename(species = Species) 

all_taxa_in_range <- bind_rows(spp_in_range, syn_to_addback, not_Actinopteri) %>%
  select(species) %>%
  unique()
# some Atlantic species end up in there because of the Arctic...
```

``` r
# add back in mammals and non-fishes
nonfish <- taxon_df %>%
  filter(class != "Actinopteri")

# keep species that are within the appropriate distribution
taxa_in_range <- taxon_df %>%
  inner_join(., all_taxa_in_range) %>%
  bind_rows(nonfish)
```

    ## Joining with `by = join_by(species)`

``` r
# blastn results that are excluded based on distributions
taxon_df %>%
  ungroup() %>%
  anti_join(., all_taxa_in_range) %>%
  select(species) %>%
  unique()
```

    ## Joining with `by = join_by(species)`

    ## # A tibble: 138 × 1
    ##    species                 
    ##    <chr>                   
    ##  1 Clupea harengus         
    ##  2 Sprattus sprattus       
    ##  3 Gadus morhua            
    ##  4 Gadus finnmarchicus     
    ##  5 Arctogadus glacialis    
    ##  6 Gadus ogac              
    ##  7 Merlangius merlangus    
    ##  8 Melanogrammus aeglefinus
    ##  9 Stichaeopsis nevelskoi  
    ## 10 Askoldia variegata      
    ## # ℹ 128 more rows

``` r
taxa_in_range 
```

    ## # A tibble: 1,384 × 9
    ## # Groups:   ASV [201]
    ##    ASV   perc_id kingdom phylum   class       order        family  genus species
    ##    <chr>   <dbl> <chr>   <chr>    <chr>       <chr>        <chr>   <chr> <chr>  
    ##  1 ASV1    100   ""      Chordata Mammalia    Carnivora    Otarii… Call… Callor…
    ##  2 ASV1     99.4 ""      Chordata Mammalia    Carnivora    Otarii… Call… Callor…
    ##  3 ASV2    100   ""      Chordata Actinopteri Clupeiformes Clupei… Clup… Clupea…
    ##  4 ASV2     99.4 ""      Chordata Actinopteri Clupeiformes Clupei… Clup… Clupea…
    ##  5 ASV3    100   ""      Chordata Actinopteri Gadiformes   Gadidae Gadus Gadus …
    ##  6 ASV3     99.4 ""      Chordata Actinopteri Gadiformes   Gadidae Gadus Gadus …
    ##  7 ASV3     99.4 ""      Chordata Actinopteri Gadiformes   Gadidae Gadus Gadus …
    ##  8 ASV3     98.8 ""      Chordata Actinopteri Gadiformes   Gadidae Gadus Gadus …
    ##  9 ASV3     98.8 ""      Chordata Actinopteri Gadiformes   Gadidae Gadus Gadus …
    ## 10 ASV3     98.2 ""      Chordata Actinopteri Gadiformes   Gadidae Bore… Boreog…
    ## # ℹ 1,374 more rows

These all look reasonable to exclude.

Now, using those data, let’s go through this procedure to identify the
appropriate taxonomic level for each sequence.

There are four categories: 1. sequences that match a single species
unambiguously (the minority)

Sequences that match multiple species are divided in three categories:
2. top matches \> 0.5% identity than second-ranked matches 3. top
matches \< 0.5% identity than second-ranked matches 4. Multiple top
matches with the same % identity

``` r
# 1. sequences that are unambiguously a single species
single_spp_seqs <- taxa_in_range %>% 
  group_by(ASV, species) %>%
  tally() %>% # how many species per ASV?
  group_by(ASV) %>%
  select(-n) %>%
  add_tally(name = "n_taxa") %>%
  filter(n_taxa == 1)
  
single_spp_seqs %>%
  left_join(., taxa_in_range) %>%
  select(perc_id, ASV, species)
```

    ## Joining with `by = join_by(ASV, species)`

    ## # A tibble: 190 × 3
    ## # Groups:   ASV [79]
    ##    perc_id ASV    species              
    ##      <dbl> <chr>  <chr>                
    ##  1   100   ASV1   Callorhinus ursinus  
    ##  2    99.4 ASV1   Callorhinus ursinus  
    ##  3   100   ASV1   Callorhinus ursinus  
    ##  4    99.4 ASV1   Callorhinus ursinus  
    ##  5    99.4 ASV101 Gadus chalcogrammus  
    ##  6    98.2 ASV101 Gadus chalcogrammus  
    ##  7    98.2 ASV101 Gadus chalcogrammus  
    ##  8    99.4 ASV106 Anoplopoma fimbria   
    ##  9    98.8 ASV106 Anoplopoma fimbria   
    ## 10   100   ASV108 Eschrichtius robustus
    ## # ℹ 180 more rows

How many species do those ASVs represent?

``` r
single_spp_seqs %>%
  ungroup() %>%
  select(species) %>%
  unique()
```

    ## # A tibble: 23 × 1
    ##    species              
    ##    <chr>                
    ##  1 Callorhinus ursinus  
    ##  2 Gadus chalcogrammus  
    ##  3 Anoplopoma fimbria   
    ##  4 Eschrichtius robustus
    ##  5 Phocoenoides dalli   
    ##  6 Clupea pallasii      
    ##  7 Trichodon trichodon  
    ##  8 Leptoclinus maculatus
    ##  9 Lycodes palearis     
    ## 10 Lycodes brevipes     
    ## # ℹ 13 more rows

``` r
# remove the single-species seqs from the dataframe and then rank the hits by % identity for the remaining seqs
seq_id_diff <- taxa_in_range %>%
  anti_join(., single_spp_seqs) %>%
  group_by(ASV, species, genus, family, order, class, phylum, kingdom) %>%
    mutate(seq_percID = max(perc_id)) %>%
    group_by(ASV, species, genus, family, order, class, phylum, kingdom, seq_percID) %>%
  summarise(max(seq_percID)) %>% # take just the top hit for each taxon (for each sequence)
  select(-`max(seq_percID)`) %>%
  ungroup() %>%
  group_by(ASV) %>%
      mutate(id_rank = rank(desc(seq_percID), ties.method = "min")) %>% # rank the taxonomic hits per sequence by % id
       mutate(top_perc = max(seq_percID)) %>% # designate the highest % id for the best taxonomic hit in each sequence (in some, but not all cases, this is 100%)   
      mutate(diff = top_perc - seq_percID) %>% # calculate the difference between the % identity of the top hit and each subsequent taxonomic hit
        arrange(diff)
```

    ## Joining with `by = join_by(ASV, species)`
    ## `summarise()` has grouped output by 'ASV', 'species', 'genus', 'family',
    ## 'order', 'class', 'phylum', 'kingdom'. You can override using the `.groups`
    ## argument.

``` r
seq_id_diff %>%
  filter(diff > 0.5)
```

    ## # A tibble: 482 × 12
    ## # Groups:   ASV [107]
    ##    ASV    species     genus family order class phylum kingdom seq_percID id_rank
    ##    <chr>  <chr>       <chr> <chr>  <chr> <chr> <chr>  <chr>        <dbl>   <int>
    ##  1 ASV164 Lethentero… Leth… Petro… Petr… Hype… Chord… ""            99.5       5
    ##  2 ASV90  Stenobrach… Sten… Mycto… Myct… Acti… Chord… ""            99.4       2
    ##  3 ASV341 Stenobrach… Sten… Mycto… Myct… Acti… Chord… ""            98.8       2
    ##  4 ASV362 Stenobrach… Sten… Mycto… Myct… Acti… Chord… ""            98.8       2
    ##  5 ASV10  Gadus macr… Gadus Gadid… Gadi… Acti… Chord… ""            99.4       2
    ##  6 ASV22  Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       2
    ##  7 ASV22  Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       2
    ##  8 ASV22  Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       2
    ##  9 ASV3   Gadus macr… Gadus Gadid… Gadi… Acti… Chord… ""            99.4       2
    ## 10 ASV39  Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       2
    ## # ℹ 472 more rows
    ## # ℹ 2 more variables: top_perc <dbl>, diff <dbl>

``` r
# salmon with 100% sequence matches?
seq_id_diff %>%
  filter(genus == "Oncorhynchus") %>%
  filter(seq_percID != 100) %>%
  arrange(ASV)
```

    ## # A tibble: 63 × 12
    ## # Groups:   ASV [15]
    ##    ASV    species     genus family order class phylum kingdom seq_percID id_rank
    ##    <chr>  <chr>       <chr> <chr>  <chr> <chr> <chr>  <chr>        <dbl>   <int>
    ##  1 ASV124 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       1
    ##  2 ASV124 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.8       2
    ##  3 ASV124 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.2       3
    ##  4 ASV124 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.2       4
    ##  5 ASV135 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       1
    ##  6 ASV135 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       1
    ##  7 ASV135 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.8       3
    ##  8 ASV135 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.2       4
    ##  9 ASV135 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.2       4
    ## 10 ASV135 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.2       4
    ## # ℹ 53 more rows
    ## # ℹ 2 more variables: top_perc <dbl>, diff <dbl>

``` r
# salmon
seq_id_diff %>%
  filter(genus == "Oncorhynchus") %>%
  arrange(ASV, id_rank) %>%
  filter(diff < 1) # remove potential hits that are >1% different
```

    ## # A tibble: 38 × 12
    ## # Groups:   ASV [15]
    ##    ASV    species     genus family order class phylum kingdom seq_percID id_rank
    ##    <chr>  <chr>       <chr> <chr>  <chr> <chr> <chr>  <chr>        <dbl>   <int>
    ##  1 ASV124 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       1
    ##  2 ASV124 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.8       2
    ##  3 ASV135 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       1
    ##  4 ASV135 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       1
    ##  5 ASV135 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.8       3
    ##  6 ASV138 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       1
    ##  7 ASV138 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       2
    ##  8 ASV138 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            98.8       3
    ##  9 ASV142 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""            99.4       1
    ## 10 ASV159 Oncorhynch… Onco… Salmo… Salm… Acti… Chord… ""           100         1
    ## # ℹ 28 more rows
    ## # ℹ 2 more variables: top_perc <dbl>, diff <dbl>

``` r
  # group_by(ASV) %>%
  # add_tally() %>%
  # filter(n == 1) %>%
  # tally() %>%
  # left_join(., seq_id_diff) %>%
  # filter(diff < 0.5)
```

``` r
seq_id_diff %>%
  filter(genus == "Gadus") %>%
  filter(seq_percID == 100) %>%
  arrange(ASV)
```

    ## # A tibble: 11 × 12
    ## # Groups:   ASV [10]
    ##    ASV    species     genus family order class phylum kingdom seq_percID id_rank
    ##    <chr>  <chr>       <chr> <chr>  <chr> <chr> <chr>  <chr>        <dbl>   <int>
    ##  1 ASV10  Gadus chal… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ##  2 ASV13  Gadus chal… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ##  3 ASV153 Gadus macr… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ##  4 ASV28  Gadus chal… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ##  5 ASV28  Gadus macr… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ##  6 ASV3   Gadus chal… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ##  7 ASV343 Gadus chal… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ##  8 ASV37  Gadus chal… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ##  9 ASV410 Gadus chal… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ## 10 ASV58  Gadus macr… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ## 11 ASV62  Gadus chal… Gadus Gadid… Gadi… Acti… Chord… ""             100       1
    ## # ℹ 2 more variables: top_perc <dbl>, diff <dbl>

Now I have the single best entry for each species for each sequence
ranked with the difference between the first and second ranked entries
calculated.

For sequences with multiple top hits, where the difference between
ranked taxa = 0, I will end up defaulting to genus- or family-level ID
(or carrying the individual species info around in some capacity). I
will do the same for any sequences where the difference betweeen the
first and second ranked taxa is \< 0.5%.

Figure out which differences are \> 0.5% and eliminate those first?

``` r
# filter out any taxa that are >1% less matching identity than the top taxonomic hit for a given sequence
to_remove_low_perc_hits <- seq_id_diff %>%
  ungroup() %>%
  group_by(ASV) %>%
  filter(diff > 0.5) # using a 0.5% threshold, there are 1719 rows in this df

keepers <- seq_id_diff %>%
  anti_join(to_remove_low_perc_hits)
```

    ## Joining with `by = join_by(ASV, species, genus, family, order, class, phylum,
    ## kingdom, seq_percID, id_rank, top_perc, diff)`

``` r
# this data frame includes only those taxonomic hits that should be considered.
# so now I need to determine whether they should be assigned to genus, family, order, etc. 

# Single-species matches now that we've filtered based on percent identity
singleton_df <- keepers %>%
  select(ASV) %>%
  tally() %>%
  filter(n == 1) %>% # these are the seqs that now have only a single match
  left_join(keepers) %>%
  select(-n) %>%
  bind_rows(single_spp_seqs) %>% # combine the single spp data
  mutate(taxonomic_level = "species") %>%
  mutate(taxon = species)
```

    ## Joining with `by = join_by(ASV)`

``` r
## Genus-level matches
# remove the singletons from the bigger df 
single_genus <- keepers %>%
  anti_join(., singleton_df, by = "ASV") %>% # at best, these should be genus-level matches
  group_by(ASV, genus) %>%
  tally() %>%
  ungroup() %>%
  group_by(ASV) %>%
  tally() %>%
  filter(n == 1) %>% # seqs that match a single genus
  select(-n) %>%
  left_join(., keepers) %>%
  mutate(taxonomic_level = "genus") %>%
  mutate(taxon = genus)
```

    ## Joining with `by = join_by(ASV)`

``` r
## Family-level matches
single_family <- keepers %>%
  anti_join(., singleton_df, by = "ASV") %>%
  anti_join(., single_genus) %>%
  group_by(ASV, family) %>%
  tally() %>%
  ungroup() %>%
  group_by(ASV) %>%
  tally() %>%
  filter(n == 1) %>% # seqs that match a single family
  select(-n) %>%
  left_join(., keepers) %>%
  mutate(taxonomic_level = "family") %>%
  mutate(taxon = family)
```

    ## Joining with `by = join_by(ASV, species, genus, family, order, class, phylum,
    ## kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(ASV)`

``` r
## Order-level matches
single_order <- keepers %>%
 anti_join(., singleton_df, by = "ASV") %>%
  anti_join(., single_genus) %>%
  anti_join(., single_family) %>%
  group_by(ASV, order) %>%
  tally() %>%
  ungroup() %>%
  group_by(ASV) %>%
  tally() %>%
  filter(n == 1) %>% # seqs that match a single order
  select(-n) %>%
  left_join(., keepers) %>%
  mutate(taxonomic_level = "order") %>%
  mutate(taxon = order)
```

    ## Joining with `by = join_by(ASV, species, genus, family, order, class, phylum,
    ## kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(ASV, species, genus, family, order, class, phylum,
    ## kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(ASV)`

``` r
## Class-level matches
single_class <- keepers %>%
  anti_join(., singleton_df, by = "ASV") %>%
  anti_join(., single_genus) %>%
  anti_join(., single_family) %>%
  anti_join(., single_order) %>%
  group_by(ASV, class) %>%
  tally() %>%
  ungroup() %>%
  group_by(ASV) %>%
  tally() %>% 
  filter(n == 1) %>% # seqs that match a single class
  select(-n) %>%
  left_join(., keepers) %>%
  mutate(taxonomic_level = "class") %>%
  mutate(taxon = class)
```

    ## Joining with `by = join_by(ASV, species, genus, family, order, class, phylum,
    ## kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(ASV, species, genus, family, order, class, phylum,
    ## kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(ASV, species, genus, family, order, class, phylum,
    ## kingdom, seq_percID, id_rank, top_perc, diff)`
    ## Joining with `by = join_by(ASV)`

``` r
# no higher level taxonomy was relevant in this dataset.
```

``` r
# combine those sub-dataframes
tax_df_w_gadidSalmon <- bind_rows(singleton_df, single_genus, single_family, single_order, single_class) %>%
  select(-id_rank, -top_perc, -diff, -seq_percID, -n_taxa) %>%
  unique()
```

Create output taxonomy data frames

``` r
uncollapsed_taxonomy <- tax_df_w_gadidSalmon %>%
  #select(-top_perc, -id_rank) %>%
  unique()

# quick look at the taxonomy before collapsing it:
uncollapsed_taxonomy %>%
  write_csv("csv_outputs/nfs_diet_mifish_inrange_uncollapsed_taxonomy.csv")

# and then collapse that down to just a single taxon per ASV
collapsed_taxonomy <- uncollapsed_taxonomy %>%
  select(ASV, taxon, taxonomic_level) %>%
  unique()

# then remove those high level taxonomic assignments that aren't informative
collapsed_taxonomy %>%
  #filter(!taxonomic_level %in% c("class", "order")) %>%
  select(-ASV) %>%
  unique()
```

    ## # A tibble: 64 × 2
    ##    taxon                    taxonomic_level
    ##    <chr>                    <chr>          
    ##  1 Gadus chalcogrammus      species        
    ##  2 Hippoglossus stenolepis  species        
    ##  3 Oncorhynchus keta        species        
    ##  4 Oncorhynchus tshawytscha species        
    ##  5 Boreogadus saida         species        
    ##  6 Blepsias bilobus         species        
    ##  7 Gadus macrocephalus      species        
    ##  8 Sebastes pinniger        species        
    ##  9 Oncorhynchus kisutch     species        
    ## 10 Lycodes brevipes         species        
    ## # ℹ 54 more rows

``` r
collapsed_taxonomy %>%
  select(-ASV) %>%
  unique() %>%
  group_by(taxonomic_level) %>%
  tally()
```

    ## # A tibble: 4 × 2
    ##   taxonomic_level     n
    ##   <chr>           <int>
    ## 1 family              5
    ## 2 genus              12
    ## 3 order               1
    ## 4 species            46

``` r
collapsed_taxonomy %>%
  write_csv("csv_outputs/nfs_diet_mifish_inrange_collapsed_taxonomy.csv")
```

Look at how many reads they got:

``` r
df <- read_csv("csv_outputs/NFS_diets_MiFish_ASVtable.csv") %>%
  rename(sample = `...1`)
```

    ## New names:
    ## Rows: 573 Columns: 430
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (1): ...1 dbl (429): ASV1, ASV2, ASV3, ASV4, ASV5, ASV6, ASV7, ASV8, ASV9,
    ## ASV10, ASV1...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
# reformat
# dealing with the positive and negative controls having a different pattern
sample_df <- df %>%
  pivot_longer(2:length(df), names_to = "ASV", values_to = "reads") %>%
  filter(reads > 0) %>%
  separate(
    sample,
    into = c("ablg", "replicate"),
    sep = "-(?=[^-]+$)",
    remove = FALSE
  )

head(sample_df)
```

    ## # A tibble: 6 × 5
    ##   sample   ablg   replicate ASV    reads
    ##   <chr>    <chr>  <chr>     <chr>  <dbl>
    ## 1 e06296-1 e06296 1         ASV1    5248
    ## 2 e06296-1 e06296 1         ASV2    8648
    ## 3 e06296-1 e06296 1         ASV7     608
    ## 4 e06296-1 e06296 1         ASV10    518
    ## 5 e06296-1 e06296 1         ASV33      3
    ## 6 e06296-1 e06296 1         ASV182    21

``` r
meta <- read_csv("../metadata/NFS_diets_metadata_20261005.csv")
```

    ## New names:
    ## Rows: 271 Columns: 44
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (14): extraction_ID, alternative_ID, source, location1, location2, locat... dbl
    ## (5): collection_year, collection_month, collection_day, preservation_vo... lgl
    ## (25): longitude, latitude, depth, sample_volume_liters, sample_collect_d...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...34`
    ## • `` -> `...35`
    ## • `` -> `...36`
    ## • `` -> `...37`
    ## • `` -> `...38`
    ## • `` -> `...39`
    ## • `` -> `...40`
    ## • `` -> `...41`
    ## • `` -> `...42`
    ## • `` -> `...43`
    ## • `` -> `...44`

``` r
# add metadata to the sample
sample_df_w_meta <- sample_df %>%
  left_join(., meta, by = c("ablg" = "extraction_ID")) %>%
  left_join(., collapsed_taxonomy, by = "ASV") %>%
  select(ablg, replicate, ASV, reads, alternative_ID, source, collection_year, location1, location2, location3, sample_type, taxon, taxonomic_level)

sample_df_w_meta %>%
  filter(sample_type == "sample") %>%
  filter(collection_year == "2025") %>%
  select(ablg, alternative_ID) %>%
  unique()
```

    ## # A tibble: 120 × 2
    ##    ablg   alternative_ID
    ##    <chr>  <chr>         
    ##  1 e06296 NFSG401ANY-9  
    ##  2 e06297 NFSG403ANY-9  
    ##  3 e06298 NFSG250ANY-9  
    ##  4 e06299 NFSG411ANY-9  
    ##  5 e06300 NFSG409ANY-9  
    ##  6 e06301 NFSG415ANY-9  
    ##  7 e06302 NFSG413ANY-9  
    ##  8 e06303 NFSG418ANY-9  
    ##  9 e06304 NFSG419ANY-9  
    ## 10 e06305 NFSG423ANY-9  
    ## # ℹ 110 more rows

# account for likely contaminants

- tag-jumping
- negative PCR controls
- field negatives

## Step 1. Account for tag-jumping by using the positive controls

subtract the proportion of reads that jumped into the positive control
samples from each environmental sample

``` r
# on a per-plate basis, what % of each ASV should be subtracted from the samples?
asv_prop_df <- sample_df %>%
  group_by(sample) %>%
  mutate(total_reads = sum(reads)) %>%
  mutate(asv_prop = reads/total_reads)
  
asv_prop_df %>%
  filter(str_detect(sample, "PC")) %>%
  ggplot(aes(x = sample, y = asv_prop, fill = ASV)) +
  geom_bar(stat = "identity") +
  theme_minimal()
```

![](02-2025-diet-analysis_files/figure-gfm/tag-jumping-by-ASVs-in-PC-1.png)<!-- -->

``` r
# based on asvs in the PC
# because there are two positive controls
# need to take the greater values of the two
asv_prop_to_remove <- asv_prop_df %>%
  filter(str_detect(sample, "PC")) %>%
  ungroup() %>%
  select(replicate, ASV, asv_prop) %>%
  group_by(replicate, ASV) %>%
  slice_max(order_by = asv_prop, n = 1) %>%
  rename(prop_to_subtract = asv_prop)

asv_prop_to_remove %>%
  arrange(replicate, ASV)
```

    ## # A tibble: 20 × 3
    ## # Groups:   replicate, ASV [20]
    ##    replicate ASV    prop_to_subtract
    ##    <chr>     <chr>             <dbl>
    ##  1 1         ASV1           0.0196  
    ##  2 1         ASV12          0.960   
    ##  3 1         ASV2           0.00158 
    ##  4 1         ASV3           0.00138 
    ##  5 1         ASV57          0.0232  
    ##  6 2         ASV1           0.0203  
    ##  7 2         ASV12          0.958   
    ##  8 2         ASV2           0.000801
    ##  9 2         ASV3           0.00184 
    ## 10 2         ASV427         0.000160
    ## 11 2         ASV57          0.0248  
    ## 12 3         ASV1           0.0208  
    ## 13 3         ASV12          0.977   
    ## 14 3         ASV2           0.000629
    ## 15 3         ASV27          0.000319
    ## 16 3         ASV3           0.00134 
    ## 17 3         ASV303         0.000839
    ## 18 3         ASV57          0.0305  
    ## 19 3         ASV7           0.000420
    ## 20 3         ASV8           0.000703

``` r
pc_filtered_sample_df <- asv_prop_df %>%
  left_join(., asv_prop_to_remove) %>%
  mutate(new_asv_prop = ifelse(!is.na(prop_to_subtract), asv_prop-prop_to_subtract, asv_prop)) %>%
  mutate(new_asv_prop = ifelse(new_asv_prop < 0, 0, new_asv_prop)) %>% # remove negative numbers
  #filter(new_asv_prop > 0) %>%
  mutate(reads_removed = prop_to_subtract*total_reads) %>%
  select(-prop_to_subtract, -asv_prop)
```

    ## Joining with `by = join_by(replicate, ASV)`

Reads are removed by proportions based on individual replicates.

Reads removed by tag-jumping:

``` r
pc_filtered_sample_df %>%
  filter(!is.na(reads_removed))
```

    ## # A tibble: 1,672 × 8
    ## # Groups:   sample [572]
    ##    sample   ablg   replicate ASV   reads total_reads new_asv_prop reads_removed
    ##    <chr>    <chr>  <chr>     <chr> <dbl>       <dbl>        <dbl>         <dbl>
    ##  1 e06296-1 e06296 1         ASV1   5248       15051      0.329          295.  
    ##  2 e06296-1 e06296 1         ASV2   8648       15051      0.573           23.8 
    ##  3 e06296-2 e06296 2         ASV1   5144       10865      0.453          221.  
    ##  4 e06296-2 e06296 2         ASV2   4997       10865      0.459            8.70
    ##  5 e06296-2 e06296 2         ASV3     60       10865      0.00368         20.0 
    ##  6 e06296-3 e06296 3         ASV1   4511        9425      0.458          196.  
    ##  7 e06296-3 e06296 3         ASV2   4238        9425      0.449            5.93
    ##  8 e06296-3 e06296 3         ASV3     80        9425      0.00715         12.6 
    ##  9 e06296-3 e06296 3         ASV7    351        9425      0.0368           3.96
    ## 10 e06297-1 e06297 1         ASV1  11016       17190      0.621          337.  
    ## # ℹ 1,662 more rows

this is a summary of the number of reads removed by ASV and sample_ID

``` r
decontaminated_1 <- pc_filtered_sample_df %>%
  dplyr::select(sample, ASV, reads_removed) %>%
  filter(!is.na(replicate) &
           !is.na(reads_removed)) %>%
  group_by(ASV) %>%
  summarise(mean_reads = mean(reads_removed),
            reads_q.05 = quantile(reads_removed, probs=0.05),
            median_q.5 = median(reads_removed),
            reads_q.95 = quantile(reads_removed, probs=0.95)) %>%
  filter(mean_reads > 0) 
```

    ## Warning: There were 573 warnings in `filter()`.
    ## The first warning was:
    ## ℹ In argument: `!is.na(replicate) & !is.na(reads_removed)`.
    ## ℹ In group 1: `sample = "NFS-01-NC-1"`.
    ## Caused by warning in `is.na()`:
    ## ! is.na() applied to non-(list or vector) of type 'closure'
    ## ℹ Run `dplyr::last_dplyr_warnings()` to see the 572 remaining warnings.

``` r
decontaminated_1  
```

    ## # A tibble: 10 × 5
    ##    ASV    mean_reads reads_q.05 median_q.5 reads_q.95
    ##    <chr>       <dbl>      <dbl>      <dbl>      <dbl>
    ##  1 ASV1       477.       214.       473.        773. 
    ##  2 ASV12    19432.     11432.     17944.      34898. 
    ##  3 ASV2        25.2        7.38      19.4        56.6
    ##  4 ASV27        6.44       3.40       5.83       11.6
    ##  5 ASV3        36.0       15.0       35.7        58.2
    ##  6 ASV303       8          8          8           8  
    ##  7 ASV427       2          2          2           2  
    ##  8 ASV57      360.       295.       352         463. 
    ##  9 ASV7         8.40       3.98       8.57       15.5
    ## 10 ASV8        17.9        9.16      16.9        29.5

## Step 2. Account for contaminants in positive and negative controls

next we will remove ASVs that only occur in controls and not in
environmental samples.

number of reads

``` r
reads_per_type_ASV <- pc_filtered_sample_df %>%
  left_join(., meta, by = c("ablg" = "extraction_ID")) %>% 
  mutate(sample_type = ifelse(is.na(sample_type) & str_detect(ablg, "NC"), "pcr_blank", sample_type)) %>%
  mutate(sample_type = ifelse(is.na(sample_type) & str_detect(ablg, "PC"), "positive_control", sample_type)) %>%
  group_by(ASV, sample_type) %>%
  summarize(TotalReadsPerASV = sum(reads, na.rm = TRUE)) %>%
  arrange(ASV)
```

    ## `summarise()` has grouped output by 'ASV'. You can override using the `.groups`
    ## argument.

what ASVs have no reads in samples, but reads in the controls?

``` r
not_in_samples <- reads_per_type_ASV %>%
  pivot_wider(names_from = "sample_type", values_from = c("TotalReadsPerASV")) %>%
  filter(is.na(sample)) %>% # here, Kim had a numeric, but mine are NA
  replace(is.na(.), 0)

not_in_samples
```

    ## # A tibble: 21 × 5
    ## # Groups:   ASV [21]
    ##    ASV    extraction_blank pcr_blank positive_control sample
    ##    <chr>             <dbl>     <dbl>            <dbl>  <dbl>
    ##  1 ASV100              317         0                0      0
    ##  2 ASV126              147         0                0      0
    ##  3 ASV131              126         0                0      0
    ##  4 ASV136              111         0                0      0
    ##  5 ASV137              105         0                0      0
    ##  6 ASV158                0        65                0      0
    ##  7 ASV235               16         0                0      0
    ##  8 ASV259               12         0                0      0
    ##  9 ASV303                0         0                8      0
    ## 10 ASV351                5         0                0      0
    ## # ℹ 11 more rows

what ASVs do have reads in samples, but more reads in the controls?

``` r
more_in_pcr_blanks <- reads_per_type_ASV %>%
  pivot_wider(names_from = "sample_type", values_from = c("TotalReadsPerASV")) %>%
  filter(sample > 1) %>%
  filter(pcr_blank > sample)
head(more_in_pcr_blanks)
```

    ## # A tibble: 0 × 5
    ## # Groups:   ASV [0]
    ## # ℹ 5 variables: ASV <chr>, extraction_blank <dbl>, pcr_blank <dbl>,
    ## #   positive_control <dbl>, sample <dbl>

``` r
more_in_extraction_blanks <- reads_per_type_ASV %>%
  pivot_wider(names_from = "sample_type", values_from = c("TotalReadsPerASV")) %>%
  filter(sample > 1) %>%
  filter(extraction_blank > sample)
head(more_in_extraction_blanks)
```

    ## # A tibble: 0 × 5
    ## # Groups:   ASV [0]
    ## # ℹ 5 variables: ASV <chr>, extraction_blank <dbl>, pcr_blank <dbl>,
    ## #   positive_control <dbl>, sample <dbl>

``` r
more_in_pc_blanks <- reads_per_type_ASV %>%
  pivot_wider(names_from = "sample_type", values_from = c("TotalReadsPerASV")) %>%
  filter(sample > 1) %>%
  filter(positive_control > sample)
head(more_in_pc_blanks)
```

    ## # A tibble: 1 × 5
    ## # Groups:   ASV [1]
    ##   ASV   extraction_blank pcr_blank positive_control sample
    ##   <chr>            <dbl>     <dbl>            <dbl>  <dbl>
    ## 1 ASV12               NA        NA            84782    267

remove these from the asv table

``` r
asv_table_filter2 <- pc_filtered_sample_df %>%
  filter(!ASV %in% not_in_samples$ASV) %>%
  select(-reads_removed, -total_reads, -new_asv_prop) %>%
 filter(!ASV %in% more_in_pcr_blanks$ASV) %>%
 filter(!ASV %in% more_in_extraction_blanks$ASV) #%>%
  #filter(!ASV %in% more_in_fb_blanks$ASV)

# Let's stick with ASVs that are just in the controls for filtering
```

## Step 3. Remove ASVs without taxonomic ID

now lets see how many of these ASVs have taxonomic IDs (these are not
final tax ids)

``` r
taxonomy <- read_csv("csv_outputs/nfs_diet_mifish_inrange_collapsed_taxonomy.csv")
```

    ## Rows: 201 Columns: 3
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (3): ASV, taxon, taxonomic_level
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
combined_df <- asv_table_filter2 %>%
  left_join(., taxonomy) 
```

    ## Joining with `by = join_by(ASV)`

``` r
combined_df %>%
  filter(is.na(taxon)) %>%
  ungroup() %>%
  select(ASV) %>%
  unique()
```

    ## # A tibble: 214 × 1
    ##    ASV   
    ##    <chr> 
    ##  1 ASV33 
    ##  2 ASV182
    ##  3 ASV262
    ##  4 ASV222
    ##  5 ASV253
    ##  6 ASV254
    ##  7 ASV263
    ##  8 ASV264
    ##  9 ASV268
    ## 10 ASV275
    ## # ℹ 204 more rows

214 without taxonomy (in range) and removed.

``` r
combined_df_noNA <- combined_df %>%
  filter(!is.na(taxon)) 

combined_df_noNA %>%
  group_by(taxon) %>%
  summarise(sum(reads))
```

    ## # A tibble: 61 × 2
    ##    taxon                      `sum(reads)`
    ##    <chr>                             <dbl>
    ##  1 Agonidae                              6
    ##  2 Ammodytes personatus              38964
    ##  3 Anoplopoma fimbria               172646
    ##  4 Antilocapra americana              1505
    ##  5 Atheresthes                         100
    ##  6 Blepsias bilobus                     77
    ##  7 Boreogadus saida                     83
    ##  8 Callorhinus ursinus             9306774
    ##  9 Clupea pallasii                 1871445
    ## 10 Cryptacanthodes aleutensis            4
    ## # ℹ 51 more rows

``` r
clean_df_filtered <- combined_df_noNA %>%
  left_join(., meta, by = c("ablg" = "extraction_ID")) %>%
  filter(sample_type == "sample" &
           collection_year == "2025") %>%
  filter(!taxon %in% c("Acipenser fulvescens", "Acipenser", "Callorhinus ursinus", "Eumetopias jubatus")) %>%
  group_by(sample, taxon) %>%
  mutate(taxon_reads = sum(reads)) %>%
  select(-ASV, -reads) %>%
  unique() %>%
  ungroup() %>%
  group_by(sample) %>%
  mutate(total_reads = sum(taxon_reads)) %>%
  select(sample, ablg, replicate, replicate, taxon_reads, taxon, total_reads, everything()) %>%
  mutate(prop_taxon = taxon_reads/total_reads)
  

clean_df_filtered %>%
  ggplot(aes(x = ablg, y = prop_taxon, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(replicate)) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90)
  ) 
```

![](02-2025-diet-analysis_files/figure-gfm/filter-controls-analyze-samples-1.png)<!-- -->

``` r
samples_to_remove <- clean_df_filtered %>%
  filter(total_reads < 1000)

clean_df_readDepth_filtered <- clean_df_filtered %>%
  anti_join(., samples_to_remove)
```

    ## Joining with `by = join_by(sample, ablg, replicate, taxon_reads, taxon,
    ## total_reads, taxonomic_level, alternative_ID, source, collection_year,
    ## collection_month, collection_day, location1, location2, location3, longitude,
    ## latitude, sample_type, depth, sample_volume_liters, sample_collect_device,
    ## filter_type, filter_size_micrometer, filter_preservation,
    ## preservation_volume_milliliters, time_of_day, project, collected_by,
    ## extraction_by, extraction_date, extraction_kit, plate_or_vial,
    ## longmire_volume_input_microliter, elution_volume_microliters,
    ## extraction_plate_ID, extraction_well_ID, dna_conc_ng_uL, additional_metadata,
    ## accompanying_data, ...34, ...35, ...36, ...37, ...38, ...39, ...40, ...41,
    ## ...42, ...43, ...44, prop_taxon)`

``` r
clean_df_readDepth_filtered %>%
  filter(taxon_reads > 10) %>%
 ggplot(aes(x = ablg, y = prop_taxon, fill = taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(replicate)) +
  theme_minimal() +
  theme(
    #legend.position = "none",
    axis.text.x = element_text(angle = 90)
  ) 
```

![](02-2025-diet-analysis_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->
Hmm. In some ways, that filtering wasn’t as effective.

Additional species to exclude: any mammals at this stage…

``` r
# what is 1500 reads of pronghorn doing in one sample?
# ringed seal?
# gray whale?
filtered_df <- clean_df_readDepth_filtered %>%
  filter(prop_taxon > 0.01) %>% # must be great than 1% of reads
  filter(!taxon %in% c("Antilocapra americana", "Pusa", "Phocoenoides dalli", "Eschrichtius robustus", "Perciformes", "Gadidae"))
```

``` r
pub_colors <- c(
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#009E73",  # green
  "#CC79A7",  # purple-pink
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#F0E442",  # yellow
  "#000000",  # black
  "#7B3294",  # deep purple
  "#A65628",  # brown
  "#1B9E77",  # teal
  "#E7298A",  # magenta
  "#66A61E",  # olive green
  "#E6AB02",  # mustard
  "#7570B3",  # slate purple
  "#E41A1C",  # red
  "#4DAF4A",  # green
  "#377EB8",  # blue
  "#984EA3",  # violet
  "#FF7F00",  # orange
  "#A6CEE3",  # light blue
  "firebrick"   # dark brown
)
```

``` r
filtered_df %>%
  ggplot(aes(x = ablg, y = prop_taxon, fill = reorder(taxon, -prop_taxon))) +
    geom_bar(stat = "identity") +
    facet_grid(rows = vars(replicate), cols = vars(location3), space = "free", scales = "free") +
    theme_bw() +
    theme(
      legend.key.size = unit(0.3, "cm"),
      axis.text.x = element_text(angle = 90, hjust = 0.9),
      axis.title.x = element_text(margin = margin(t = 10))
    )  +
  labs(
    fill = "taxon",
    title = "2025 St. Paul NFS diets",
    x = "sample", 
    y = "proportion of reads"
  ) +
  scale_fill_manual(values = pub_colors)
```

![](02-2025-diet-analysis_files/figure-gfm/plot-using-rookery-info-1.png)<!-- -->

``` r
ggsave("pdf_outputs/2025_stpaul_filtered_readProp.png", width = 10, height = 4)
```

Take a closer look at salmon:

``` r
filtered_df %>%
  filter(str_detect(taxon, "Oncorhynchus") | str_detect(taxon, "Gadus")) %>%
  ggplot(aes(x = alternative_ID, y = prop_taxon, fill = reorder(taxon, -prop_taxon))) +
    geom_bar(stat = "identity") +
    facet_grid(rows = vars(replicate), cols = vars(location3), space = "free", scales = "free") +
    theme_bw() +
    theme(
           legend.key.size = unit(0.3, "cm"),
      axis.text.x = element_text(angle = 90, hjust = 0.9),
      axis.title.x = element_text(margin = margin(t = 10))
    )  +
  labs(
    fill = "taxon",
    title = "Salmon in 2025 St. Paul NFS diets",
    x = "sample", 
    y = "proportion of reads"
  )  +
  scale_fill_manual(values = pub_colors)
```

![](02-2025-diet-analysis_files/figure-gfm/salmon-in-diets-1.png)<!-- -->

``` r
ggsave("pdf_outputs/2025_stpaul_filtered_salmon_gadids.png", width = 8, height = 4)
```

``` r
filtered_df %>%
  ggplot(aes(x = ablg, y = prop_taxon, fill = reorder(taxon, -prop_taxon))) +
    geom_bar(stat = "identity") +
    facet_grid(rows = vars(replicate), cols = vars(source), space = "free", scales = "free") +
    theme_bw() +
    theme(
      legend.key.size = unit(0.3, "cm"),
      axis.text.x = element_text(angle = 90, hjust = 0.9),
      axis.title.x = element_text(margin = margin(t = 10))
    )  +
  labs(
    fill = "taxon",
    title = "2025 St. Paul NFS diets, fecal loops vs. scats",
    x = "sample", 
    y = "proportion of reads"
  ) +
  scale_fill_manual(values = pub_colors)
```

![](02-2025-diet-analysis_files/figure-gfm/look-at-fecal-loops-vs-scats-1.png)<!-- -->

``` r
ggsave("pdf_outputs/2025_stpaul_filtered_readProp_source.png", width = 10, height = 4)
```

``` r
# filtered_df %>%
#   ungroup() %>%
#   select(taxon) %>%
#   unique() %>%
#   write_csv("csv_outputs/taxa_names.csv")

# read in data with common names
common_names <- read_csv("csv_outputs/taxa_w_common_names.csv")
```

    ## Rows: 22 Columns: 2
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (2): taxon, common_name
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# join those to the filtered df
```

``` r
# deal with the duplicate samples from scats
filtered_df %>%
  ungroup() %>%
  select(ablg, alternative_ID, taxon, source) %>%
  unique() %>%
  mutate(
    sample_id = if_else(
      str_detect(alternative_ID, "^NFS"),
      str_remove(alternative_ID, "[A-Za-z]$"),
      alternative_ID
    )
  ) %>%
  select(-alternative_ID, -ablg) %>% # assuming the sample_id is unique to each physical sample.
  unique() %>%
  group_by(taxon, source) %>%
  tally(name = "n_samples_per_taxon") %>%
  arrange(desc(n_samples_per_taxon)) %>%
  write_csv("csv_outputs/2025_prib_filtered_taxon_FO.csv")
```

Now looking at the rookeries according to the English Bay/East
categories in the ESR.

Katie Luxa told me that Vostochni and Polovina Cliffs –\> SP East
Zapadni Reef –\> SP English Bay

So given that… Exclude the fecal loops for the time being. Compare the
scats across the two areas.

``` r
# this version of calculating FO doesn't account for the duplicate sample IDs across rookeries...
n_sample_by_region <- filtered_df %>%
  filter(source != "fecal loop") %>%
  mutate(region = ifelse(location3 == "Zapadni Reef", "SP English Bay", "SP East")) %>%
  ungroup() %>%
  select(ablg, alternative_ID, taxon, region) %>%
  unique() %>%
  mutate(
    sample_id = if_else(
      str_detect(alternative_ID, "^NFS"),
      str_remove(alternative_ID, "[A-Za-z]$"),
      alternative_ID
    )
  ) %>%
  select(-alternative_ID, -ablg) %>% # assuming the sample_id is unique to each physical sample.
  unique() %>%
  group_by(sample_id, region) %>%
  tally(name = "n_species_per_sample") %>%
  group_by(region) %>%
  add_tally(name = "n_samples_by_region") %>%
  select(region, n_samples_by_region) %>%
  ungroup()

# add that info onto the earlier dataframe that includes species
fo_data <- filtered_df %>%
  filter(source != "fecal loop") %>%
  mutate(region = ifelse(location3 == "Zapadni Reef", "SP English Bay", "SP East")) %>%
  ungroup() %>%
  select(ablg, alternative_ID, taxon, region) %>%
  unique() %>%
  mutate(
    sample_id = if_else(
      str_detect(alternative_ID, "^NFS"),
      str_remove(alternative_ID, "[A-Za-z]$"),
      alternative_ID
    )
  ) %>%
  select(-alternative_ID, -ablg) %>% # assuming the sample_id is unique to each physical sample.
  unique() %>%
  group_by(taxon, region) %>%
  tally(name = "n_samples_per_taxon") %>%
  arrange(desc(n_samples_per_taxon)) %>%
  left_join(., n_sample_by_region) %>%
  unique() %>%
  mutate(freq = n_samples_per_taxon/n_samples_by_region*100)
```

    ## Joining with `by = join_by(region)`

    ## Warning in left_join(., n_sample_by_region): Detected an unexpected many-to-many relationship between `x` and `y`.
    ## ℹ Row 1 of `x` matches multiple rows in `y`.
    ## ℹ Row 2 of `y` matches multiple rows in `x`.
    ## ℹ If a many-to-many relationship is expected, set `relationship =
    ##   "many-to-many"` to silence this warning.

But, if we wanted to take into account the relative proportion of the
diet items, how might that change things?

## hard parts data

To compare sample-to-sample

``` r
hp <- read_csv("../data/hard_parts/NFS Diet 2025.csv")
```

    ## Rows: 150 Columns: 12
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: ","
    ## chr (10): species, sitename, complex, samplelabel, type, dna, keep, firstgro...
    ## dbl  (2): date, speciesid
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
# lots of formatting required.
hp_w_dna <- hp %>%
  mutate(common_name = case_when(
    primary == "Gadus" ~ "unid. gadid",
    primary == "Herring" ~ "Pacific herring",
    primary == "Pollock" ~ "Walleye pollock",
    primary == "Salmon" ~ "Salmon sp.",
    primary == "Sandfish" ~ "Pacific sandfish",
    firstgroupname == "Pleurogrammus monopterygius/Anoplopoma fimbria" ~ "Atka mackerel or sablefish",
    primary == "Sandlance" ~ "Pacific sandlance",
    firstgroupname == "Anoplopoma fimbria" ~ "Sablefish",
    firstgroupname == "Hexagrammos sp." ~ "Greenling sp.",
    .default = primary
    
  )
           ) %>%
  filter(dna != "no") %>%
  #select(sitename, dna) %>% # the sitename plus dna id is unique.
  #unique() %>%
  separate(sitename, into = c("loc1", "location3"), 
           sep = " \\(",
    remove = TRUE) %>%
    mutate(location3 = str_remove(location3, "\\)$"),
    across(c(loc1, location3), str_to_title)) %>%
  mutate(scat_id = paste0(location3,"_",dna)) # create a unique identifier
  
tmp1 <- hp_w_dna %>%
  select(scat_id, keep) %>%
  unique() %>%
  mutate(analysis = "hard_parts")
```

Why does this have 51 samples? Looks like a couple of fecal loops in
there…

``` r
# explicitly deal with the scat subsamples
mb_subsamples <- filtered_df %>%
  ungroup() %>%
  filter(source != "fecal loop") %>%
  separate(alternative_ID, into = c("dna", "subsample"),
            sep = "(?<=\\d)(?=[A-Za-z])") %>%
  select(location3, dna, subsample, taxon, prop_taxon) %>%
  unique() %>%
  group_by(location3, dna, subsample, taxon) %>%
  mutate(mean_prop = mean(prop_taxon)) %>%
  select(-prop_taxon) %>%
  unique() %>%
  ungroup() %>%
  mutate(scat_id = paste0(location3,"_",dna)) %>%
  mutate(region = ifelse(location3 == "Zapadni Reef", "SP English Bay", "SP East"))


mb_subsamples %>%
  ggplot(aes(x = scat_id, y = mean_prop, fill = reorder(taxon, -mean_prop))) +
  geom_bar(stat = "identity") +
  facet_grid(rows = vars(subsample), cols = vars(region), space = "free", scales = "free", labeller = label_both) +
  theme_bw() +
    theme(
      legend.key.size = unit(0.4, "cm"),
      #axis.text.x = element_text(angle = 90, hjust = 0.9),
      axis.text.x = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10)),
      legend.position = "bottom"
    )  +
  labs(
    fill = "",
    title = "2025 St. Paul Northern Fur Seal diets\nDNA metabarcoding of scats",
    x = "sample", 
    y = "proportion of reads"
  ) +
  scale_fill_manual(values = pub_colors)
```

![](02-2025-diet-analysis_files/figure-gfm/getting-data-oriented-for-ESR-1.png)<!-- -->

``` r
# better to take the mean across all replicates and subsamples at the same time!
uniq_scat_sample_df <- filtered_df %>%
  ungroup() %>%
  filter(source != "fecal loop") %>%
  separate(alternative_ID, into = c("dna", "subsample"),
            sep = "(?<=\\d)(?=[A-Za-z])") %>%
  select(location3, dna, taxon, prop_taxon) %>%
  unique() %>%
  group_by(location3, dna, taxon) %>%
  mutate(mean_prop = mean(prop_taxon)) %>%
  select(-prop_taxon) %>%
  unique() %>%
  ungroup() %>%
  filter(mean_prop > 0.02) %>% # greater than 2% of reads
  # and then rescale
  group_by(location3, dna) %>%
  mutate(tot_prop = sum(mean_prop)) %>%
  mutate(taxon_prop = mean_prop/tot_prop) %>%
  mutate(scat_id = paste0(location3,"_",dna)) %>%
  mutate(region = ifelse(location3 == "Zapadni Reef", "SP English Bay", "SP East")) %>%
  unique() %>%
  ungroup() %>%
  left_join(., common_names) 
```

    ## Joining with `by = join_by(taxon)`

``` r
uniq_scat_sample_df %>%
    mutate(sp_grp = case_when(
    common_name == "Pacific herring" ~ 1,
    common_name == "Walleye pollock" ~ 2,
    common_name == "Prowfish" ~ 4,
    common_name == "Salmon sp." ~ 5,
    common_name == "Chum salmon" ~ 5,
    common_name == "Chinook salmon" ~ 5,
    common_name == "Pink salmon" ~ 5,
    common_name == "Sockeye salmon" ~ 5,
    common_name == "Kelp greenling" ~ 6,
    common_name == "Greenling sp." ~ 6,
    common_name == "Sablefish" ~ 7,
    .default = 8
    
  )
  ) %>%
  ggplot(aes(x = scat_id, y = taxon_prop, fill = reorder(common_name, sp_grp))) +
  geom_bar(stat = "identity") +
  facet_grid(cols = vars(region), space = "free", scales = "free") +
  theme_bw() +
    theme(
      legend.key.size = unit(0.7, "cm"),
      #axis.text.x = element_text(angle = 90, hjust = 0.9),
      axis.text.x = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10)),
      legend.position = "bottom"
    )  +
  labs(
    fill = "",
    title = "2025 St. Paul Northern Fur Seal diets\nDNA metabarcoding of scats",
    x = "sample", 
    y = "proportion of sequencing reads"
  ) +
  scale_fill_manual(values = pub_colors) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(nrow = 4))
```

![](02-2025-diet-analysis_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

``` r
ggsave("pdf_outputs/2025_stpaul_filtered_readProp_complex.png", width = 8, height = 5)
```

``` r
n_scats_per_region <- uniq_scat_sample_df %>%
  select(region, scat_id) %>%
  unique() %>%
  group_by(region) %>%
  tally(name = "n_scats_region")
```

35 total scats (excluding fecal loops)

``` r
dna_fo_df <- uniq_scat_sample_df %>%
  group_by(region, common_name) %>%
  tally() %>%
  left_join(., n_scats_per_region) %>%
  mutate(freq = n/n_scats_region, type = "DNA") 
```

    ## Joining with `by = join_by(region)`

``` r
dna_fo_df %>%
  ggplot(aes(x = freq, y = reorder(common_name, freq), fill = region)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  theme_minimal() +
  labs(
    x = "Frequency of occurrence (%)",
    y = "",
    title = "St. Paul - 2025 metabarcoding diets"
  ) +
  theme(
    axis.text.y = element_text(vjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 5))
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("navyblue", "goldenrod1"))
```

![](02-2025-diet-analysis_files/figure-gfm/look-at-fo-with-that-df-1.png)<!-- -->

``` r
ggsave("pdf_outputs/2025_stpaul_filtered_FO.png", width = 5.5, height = 4)
```

``` r
# comparing samples for hard parts that say they were analyzed for DNA with our mb data
uniq_scat_sample_df %>%
  select(region, scat_id) %>%
  unique() %>%
  full_join(., tmp1) %>%
  # ignore fecal loop samples because those shouldn't be included ad-hoc
  filter(!str_detect(scat_id, "ANY")) %>%
  mutate(dna_analysis = ifelse(is.na(region), "exclude", "keep")) %>%
  rename(hp_analysis = keep) %>%
  group_by(hp_analysis, dna_analysis) %>%
  tally()
```

    ## Joining with `by = join_by(scat_id)`

    ## # A tibble: 4 × 3
    ## # Groups:   hp_analysis [2]
    ##   hp_analysis dna_analysis     n
    ##   <chr>       <chr>        <int>
    ## 1 Exclude     exclude         13
    ## 2 Exclude     keep             5
    ## 3 Keep        exclude          1
    ## 4 Keep        keep            32

Of the samples that were collected for DNA analysis, 14 did not have the
minimum data quality required for analysis, 13 of which were also
excluded from hard parts results.

Now try combining hard parts with metabarcoding data…

``` r
n_hp_scats_complex <- hp_w_dna %>%
  filter(keep != "Exclude") %>%
  select(complex, scat_id) %>%
  unique() %>%
  group_by(complex) %>%
  tally(name = "n_scats_complex")
  
  
hp_fo_df <- hp_w_dna %>%
  filter(keep != "Exclude") %>%
  select(scat_id, complex, common_name, firstgroupname, secondgroupname) %>%
  group_by(complex, common_name) %>%
  tally() %>%
  left_join(., n_hp_scats_complex) %>%
  mutate(freq = n/n_scats_complex, type = "hard_parts") %>%
  rename(region = complex) %>%
  mutate(region = ifelse(region == "East", "SP East", "SP English Bay"))
```

    ## Joining with `by = join_by(complex)`

``` r
# combine dna and hard parts
tmp2 <- dna_fo_df %>%
  bind_rows(hp_fo_df) %>%
  filter(common_name != "Polychaeta") %>%
  mutate(sp_grp = case_when(
    common_name == "Pacific herring" ~ 1,
    common_name == "Walleye pollock" ~ 2,
    common_name == "unid. gadid" ~ 3,
    common_name == "Prowfish" ~ 4,
    common_name == "Salmon sp." ~ 5,
    common_name == "Chum salmon" ~ 5,
    common_name == "Chinook salmon" ~ 5,
    common_name == "Pink salmon" ~ 5,
    common_name == "Sockeye salmon" ~ 5,
    common_name == "Kelp greenling" ~ 6,
    common_name == "Greenling sp." ~ 6,
    common_name == "Sablefish" ~ 7,
    common_name == "Atka mackerel or sablefish" ~ 7,
    .default = 8
    
  )
  )

tmp2 %>%
  ggplot(aes(x = freq, y = reorder(common_name, -sp_grp), fill = type)) +
  geom_bar(stat = "identity", position = position_dodge2(preserve = "single")) +
  facet_grid(cols = vars(region), scales = "free", space = "free") +
  theme_bw() +
  labs(
    x = "Frequency of occurrence (%)",
    y = "",
    title = "Northern Fur Seal diet analysis:\nhard parts & DNA metabarcoding",
    fill = "method"
  ) +
  theme(
    axis.text.y = element_text(vjust = 0.5),
    axis.title.x = element_text(margin = margin(t = 5))
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_fill_manual(values = c("cornflowerblue","navyblue"))
```

![](02-2025-diet-analysis_files/figure-gfm/first-deal-with-hard-parts-fo-1.png)<!-- -->

``` r
ggsave("pdf_outputs/2025_stpaul_filtered_FO_comparison.png", width = 6.5, height = 4.5)
```

Now the challenge is that Tonya wants families grouped together… which
probably needs to be done manually.
