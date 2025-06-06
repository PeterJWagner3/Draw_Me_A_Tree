#### Setup Program ####
tree_drawing_folder <- "~/Documents/R_Projects/Tree_Drawing/";
setwd(tree_drawing_folder);
common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";	# directory to folder where you keep common source
franky <- "Franklin Gothic Medium";
franky_italic <- "Franklin Gothic Medium Italic";

#install.packages("png")
library(png);
library(readxl);
library(zoo);
source(paste(common_source_folder,"Chronos.r",sep="")); 							
source(paste(common_source_folder,"Data_Downloading_v4.r",sep="")); 				# Stuck? Use source(file.choose()) & grab it!
source(paste(common_source_folder,"Disparity.r",sep="")); 				# Stuck? Use source(file.choose()) & grab it!
source(paste(common_source_folder,"Geography.r",sep="")); 							
source(paste(common_source_folder,"General_Plot_Templates.r",sep=""));
source(paste(common_source_folder,"Nexus_File_Routines.r",sep=""));
source(paste(common_source_folder,"Occurrence_Data_Routines.r",sep=""));						
source(paste(common_source_folder,"paleophylogeny_routines.r",sep=""));						
source(paste(common_source_folder,"Stratigraphy.r",sep=""));						
source(paste(common_source_folder,"Wagner_Stats_and_Probability_101.r",sep=""));	
source(paste(common_source_folder,"Wagner_Kluges.r",sep="")); 						

load(paste(data_for_R_folder,"Gradstein_2020_Augmented.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
load(paste(data_for_R_folder,"Rock_Unit_Database.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
options(warn=0);
load(paste(data_for_R_folder,"Paleobiology_Database.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
load(paste(data_for_R_folder,"PaleoDB_Edits.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
options(warn=1);

pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
pbdb_taxonomy <- rbind(pbdb_taxonomy,paleodb_fixes$paleodb_taxonomy_edits[!paleodb_fixes$paleodb_taxonomy_edits$taxon_no %in% pbdb_taxonomy$taxon_no,]);
time_scale <- gradstein_2020_emended$time_scale;

gap <- INAP <- -22;
UNKNOWN <- missing <- -11;
polymorphs <- TRUE	# if false, then these are converted to unknowns

# set parameters & files for printing tree ####
ingroup <- "Acanthoceratidae";
mctree_file <- paste(tree_drawing_folder,"post_rj_Acanthoceratidae_250K_MCC.tre",sep="");
range_file <- paste(tree_drawing_folder,"Acanthoceratidae_Range_Data.xlsx",sep="");
taxon_information_file <- "";
split_outgroups <- FALSE;
ingroup_rank <- pbdb_taxonomy$accepted_rank[match(ingroup,pbdb_taxonomy$taxon_name)]
#rank_to_plot <- "subfamily";
strat_rank <- "Substage";		# I recommend either "Stage" or "Substage"
stratigraphic_scale <- "International";		# I recommend either "International" or "Stage Slice"
study_end <- "Turonian";
study_ma_ub <- time_scale$ma_ub[time_scale$interval==study_end];

# load stratigraphic data & set it up for figures ####
if (gsub(".csv","",range_file)!=range_file) {
	fuzzy_ranges <- read.csv(range_file,header=TRUE,stringsAsFactors = FALSE);
	} else if (gsub(".xlsx","",range_file)!=range_file)  {
	fuzzy_ranges <- as.data.frame(read_xlsx(range_file));
  }
colnames(fuzzy_ranges)[colnames(fuzzy_ranges)=="species"] <- "taxon";
colnames(fuzzy_ranges) <- tolower(colnames(fuzzy_ranges));
fuzzy_ranges$taxon <- gsub("_"," ",fuzzy_ranges$taxon);
fuzzy_ranges$la_ub[fuzzy_ranges$la_ub<study_ma_ub] <- fuzzy_ranges$la_lb[fuzzy_ranges$la_lb<study_ma_ub] <- study_ma_ub;
notu <- nrow(fuzzy_ranges);

if (taxon_information_file=="")  {
	total_ingroup <- accersi_daughter_taxa(parent_taxon=ingroup,pbdb_taxonomy);
	if (ingroup_rank %in% colnames(pbdb_taxonomy))	{
		ingroup_taxonomy <- pbdb_taxonomy[pbdb_taxonomy[,ingroup_rank] %in% ingroup,];
		} else {
		ingroup_taxonomy <- pbdb_taxonomy[pbdb_taxonomy$taxon_name==ingroup,];
		ingroup_taxonomy <- rbind(ingroup_taxonomy,pbdb_taxonomy[pbdb_taxonomy$parent_no %in% ingroup_taxonomy$accepted_no[1],]);
		for (tx in 1:notu)
			ingroup_taxonomy <- rbind(ingroup_taxonomy,pbdb_taxonomy[pbdb_taxonomy$taxon_name==fuzzy_ranges$taxon[tx],])
		ingroup_taxonomy <- add_higher_taxon_to_pbdb_taxonomy(taxon_rank=ingroup_rank,paleodb_taxonomy = ingroup_taxonomy);
		}
	ingroup_taxonomy <- pretend_subgenera_are_genera(pbdb_taxonomic_data=ingroup_taxonomy);
	ingroup_taxonomy <- accersi_faux_genus_species_combos_for_species_in_subgenera(ingroup_taxonomy);
	otu_information <- data.frame(taxon=fuzzy_ranges$taxon);
	otu_information$genus <- ingroup_taxonomy$genus[match(otu_information$taxon,ingroup_taxonomy$taxon_name)];
	otu_information$family <- ingroup_taxonomy$family[match(otu_information$taxon,ingroup_taxonomy$taxon_name)];
  otu_information$order <- ingroup_taxonomy$order[match(otu_information$taxon,ingroup_taxonomy$taxon_name)];
  otu_information$class <- ingroup_taxonomy$class[match(otu_information$taxon,ingroup_taxonomy$taxon_name)];
  otu_information$taxon_type <- rep("outgroup",notu);
  dt <- 0;
  while (dt < nrow(total_ingroup))  {
    dt <- dt+1;
    otu_information$taxon_type[which(otu_information==total_ingroup$daughter[dt],arr.ind=TRUE)[,1]] <- "ingroup";
    }
  } else if (gsub(".csv","",taxon_information_file)!=taxon_information_file) {
	otu_information <- read.csv(taxon_information_file,header=TRUE,stringsAsFactors = FALSE);
	} else if (gsub(".xlsx","",taxon_information_file)!=taxon_information_file) {
	otu_information <- as.data.frame(read_xlsx(taxon_information_file));
	}

otu_taxonomy <- ingroup_taxonomy[ingroup_taxonomy$taxon_name %in% otu_information$taxon,];

for (nn in 1:notu)	{
	otu_sites <- pbdb_sites[pbdb_sites$collection_no %in% pbdb_finds$collection_no[pbdb_finds$accepted_no %in% otu_taxonomy$accepted_no[nn]],];
	}

#ingroup_taxonomy[gsub("niciensis","",ingroup_taxonomy$taxon_name)!=ingroup_taxonomy$taxon_name,]
#if (sum(is.na(otu_information$genus))>0)	{
#	nahs <- (1:notu)[is.na(otu_information$genus)]
#	ptaxa <- nrow(pbdb_taxonomy);	
#	for (nn in 1:length(nahs))	{
#		this_genus <- divido_genus_names_from_species_names(otu_information$taxon[nahs[nn]]);
#		if (is.subgenus(this_genus))
#			this_genus <- divido_subgenus_names_from_genus_names(this_genus);
#		this_species <- divido_species_epithets(otu_information$taxon[nahs[nn]]);
#		poss_gen <- (1:ptaxa)[gsub(this_genus,"",pbdb_taxonomy$taxon_name)!=pbdb_taxonomy$taxon_name];
#		poss_spc <- (1:ptaxa)[gsub(this_species,"",pbdb_taxonomy$taxon_name)!=pbdb_taxonomy$taxon_name];
#		otu_information$genus[nahs[nn]] <- pbdb_taxonomy[poss_spc[poss_spc %in% poss_gen],]$genus[1];
#		otu_information$family[nahs[nn]] <- pbdb_taxonomy[poss_spc[poss_spc %in% poss_gen],]$family[1];
#		otu_information$order[nahs[nn]] <- pbdb_taxonomy[poss_spc[poss_spc %in% poss_gen],]$order[1];
#		otu_information$class[nahs[nn]] <- pbdb_taxonomy[poss_spc[poss_spc %in% poss_gen],]$class[1];
#		if (ingroup %in% otu_information[nahs[nn],])
#			otu_information$taxon_type[nahs[nn]] <- "ingroup";
#		}
#	}

#colnames(otu_information)[colnames(otu_information)=="species"] <- "taxon";
#otu_information$genus <- pbdb_taxonomy$genus[match(otu_information$taxon,pbdb_taxonomy$taxon_name)];
#otu_information$family <- pbdb_taxonomy$family[match(fuzzy_ranges$taxon,pbdb_taxonomy$taxon_name)];
#otu_information$order <- pbdb_taxonomy$order[match(fuzzy_ranges$taxon,pbdb_taxonomy$taxon_name)];
#otu_information$class <- pbdb_taxonomy$class[match(fuzzy_ranges$taxon,pbdb_taxonomy$taxon_name)];
if (is.null(otu_information$taxon_type))  {
	otu_information$taxon_type <- rep("outgroup",notu);
	if (length(ingroup_rank)==1)	{
		outgr_col <- match(ingroup_rank,colnames(otu_information));
		otu_information$taxon_type[otu_information[,outgr_col]==ingroup] <- "ingroup";
		} else if (length(ingroup_rank)==0)	{
		ingroup_info <- pbdb_taxonomy[pbdb_taxonomy$taxon_name %in% ingroup & pbdb_taxonomy$flags=="",];
		otu_taxonomy <- rbind(pbdb_taxonomy[pbdb_taxonomy$taxon_name %in% otu_information$taxon,],ingroup_info)
		otu_taxonomy <- add_higher_taxon_to_pbdb_taxonomy(taxon_rank=ingroup_taxonomy$accepted_rank,otu_taxonomy)
#		pbdb_taxonomy[pbdb_taxonomy$parent_name=="Eublastoidea",]
		}
	}

finest_chronostrat <- time_scale[time_scale$scale==stratigraphic_scale,];
finest_chronostrat <- finest_chronostrat[finest_chronostrat$chronostratigraphic_rank==strat_rank & finest_chronostrat$interval_sr=="",];
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];
finest_chronostrat <- finest_chronostrat[finest_chronostrat$ma_lb>=min(fuzzy_ranges$la_ub),];

# get environmental data 
for (sp in 1:notu)	{
	
	ingroup_taxonomy$accepted_no[ingroup_taxonomy$taxon_name %in% otu_names[sp]]
	}

#### Get Tree Information for Plotting ####
mctree <- scan(file=mctree_file,what=character(),sep="\n");
mctree <- gsub("\t","",mctree)
begin_taxa <- match("begin taxa;",tolower(mctree));
notu <- as.numeric(strsplit(gsub(";","",mctree[begin_taxa+1]),split = "=")[[1]][2]);

taxlabels <- match("taxlabels",tolower(mctree));
otu_names_nex <- mctree[(taxlabels+1):length(mctree)];
otu_names_nex <- otu_names_nex[1:(match(";",otu_names_nex)-1)];
otu_names <- gsub("_"," ",otu_names_nex);
# remove any taxa from ranges that are not in the tree;
fuzzy_ranges <- fuzzy_ranges[fuzzy_ranges$taxon %in% otu_names,];
otu_information <- otu_information[otu_information$taxon %in% otu_names,];
# make sure that ranges & taxa are in the same order!
fuzzy_ranges <- fuzzy_ranges[match(otu_names,fuzzy_ranges$taxon),];
otu_information <- otu_information[match(otu_names,otu_information$taxon),];
notu <- nrow(fuzzy_ranges);
# use this for anagenesis.
min_ranges <- data.frame(taxon=fuzzy_ranges$taxon,fa_mr=rep(0,notu),la_mr=rep(0,notu));
for (sp in 1:notu)  {
  if (abs(fuzzy_ranges$fa_ub[sp])>=abs(fuzzy_ranges$la_lb[sp]))  {
    min_ranges$fa_mr[sp] <- fuzzy_ranges$fa_ub[sp];
    min_ranges$la_mr[sp] <- fuzzy_ranges$la_lb[sp];
    }
  }

tree_line <- 1+match("begin trees;",tolower(mctree));
if (tolower(mctree[tree_line])=="translate")	{
	tree_line <- tree_line+notu+2;
	beast <- TRUE;
	} else	BEAST <- FALSE;
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&R\\] ","",mctree[tree_line]);
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&R\\]","",mctree[tree_line]);
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&U\\]","",newick_string_full);
newick_string_full <- tree_info <- gsub(" ","",newick_string_full);
if (gsub(" = \\[&","",newick_string_full)!=newick_string_full)	{
	newick_atomized <- strsplit(newick_string_full,"")[[1]];
	st <- match("(",newick_atomized);
	newick_string_full <- tree_info <- paste(newick_atomized[st:length(newick_atomized)],collapse="");
	}

if (!beast) {
	tree_output <- read_newick_string_mcmc(newick_string_full,otu_names);
	} else	{
	tree_output <- read_newick_string_mcmc_beast(newick_string_full,otu_names);
	}
mctree_file2 <- paste(tree_drawing_folder,"Strophomenoid_Punctuated_MCT.tre",sep="");
mctree2 <- scan(file=mctree_file2,what=character(),sep="\n")[1];
mctree2_info <- read_newick_string_with_branch_lengths(newick_string=mctree2);
clade_posteriors <- tree_output$clade_posteriors;
tree_output$ancestral <- tree_output$ancestral[!is.na(names(tree_output$ancestral))];
sampled_ancestors <- tree_output$ancestral;

notu <- match(-1,tree_output$vector_tree)-1;
branch_durations2 <- mctree2_info$branch_lengths;
vector_tree2 <- mctree2_info$vector_tree;
htu_heights2 <- mctree2_info$htu_heights;
nhtu <- length(vector_tree2);

br_lngth_history <- array(0,dim=c(notu,1+max(htu_heights2)));
rownames(br_lngth_history) <- names(vector_tree2)[1:notu];
for (sp in 1:notu)	{
	anc_node <- vector_tree2[sp];
	pd <- 1;
	br_lngth_history[sp,pd] <- branch_durations2[sp];
#	pd <- pd+1;
#	br_lngth_history[sp][pd] <- branch_durations[anc_node];
	while (anc_node!=-1)	{
		pd <- pd+1;
		br_lngth_history[sp,pd] <- branch_durations2[anc_node];
		anc_node <- vector_tree2[anc_node];
		}
	}
br_lngth_history <- br_lngth_history[,colSums(br_lngth_history)>0];

faux_recent <- min(fuzzy_ranges$fa_ub);
oldest_fa <- -(max(fuzzy_ranges$fa_ub)+max(fuzzy_ranges$fa_lb))/2;
rescaled_fa_bounds <- fuzzy_ranges;
rescaled_fa_bounds$la_lb <- rescaled_fa_bounds$la_ub <- NULL;
rescaled_fa_bounds$fa_lb <- rescaled_fa_bounds$fa_lb - faux_recent;
rescaled_fa_bounds$fa_ub <- rescaled_fa_bounds$fa_ub - faux_recent;

basal_taxon <- fuzzy_ranges$taxon[match(max(fuzzy_ranges$fa_lb),fuzzy_ranges$fa_lb)];
basal_date_raw <- max(((tree_output$hpd_age$lb+tree_output$hpd_age$ub)/2)[1:notu]);	# base from Bayesian tree before rescaling 0 to actual upper bound
basal_taxon_no <- match(basal_date_raw,((tree_output$hpd_age$lb+tree_output$hpd_age$ub)/2)[1:notu]);
#date_rescale <- ((fuzzy_ranges$fa_lb[match(basal_taxon,fuzzy_ranges$taxon)]+fuzzy_ranges$fa_ub[match(basal_taxon,fuzzy_ranges$taxon)])/2)-basal_date_raw;
#faux_recent <- min(fuzzy_ranges$fa_ub);
hpd_age <- -(tree_output$hpd_age+faux_recent);
hpd_age$md <- (hpd_age$lb+hpd_age$ub)/2;
#if (max(hpd_age$md[1:notu] - -abs(fuzzy_ranges$fa_lb))>0)	{
#	adj <- max(hpd_age$md[1:notu] - -abs(fuzzy_ranges$fa_lb));
#	hpd_age <- hpd_age+adj;
#	}
#paracrinoid_hpd_age <- hpd_age;
clade_onset <- max(abs(hpd_age));
finest_chronostrat <- finest_chronostrat[finest_chronostrat$ma_ub<clade_onset,];

#hpd_age <- data.frame(-min(fuzzy_ranges$fa_lb)-tree_output$hpd_age);
branch_durations <- abs(tree_output$branch_durations);

#gsub("_"," ",names(branch_durations)[1:notu])[!gsub("_"," ",names(branch_durations)[1:notu]) %in% fuzzy_ranges$taxon]

#branch_durations[(1:notu)[sampled_ancestors==1]] <- 0;
vector_tree <- tree_output$vector_tree;
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree=tree_output$vector_tree);
nNodes <- nrow(mat_tree);
lowest_taxon <- min(which(mat_tree %in% (1:notu),arr.ind = TRUE)); # get the OTU that is closest to the base of the tree
outgp <- mat_tree[lowest_taxon,mat_tree[lowest_taxon,]<=notu];
#-max(abs(hpd_age$md[outgp]))

#### Plot like a Lannister ####
if (variable_to_plot=="taxonomy")  {
  plot_col <- match(tolower(rank_to_plot),tolower(colnames(otu_information)));
  } else if (variable_name=="geography")  {
	fuzzy_ranges$region[fuzzy_ranges$region==""] <- fuzzy_ranges$realm[fuzzy_ranges$region==""];
	devonian_regions <- unique(fuzzy_ranges$region[fuzzy_ranges$fa_lb <= time_scale$ma_lb[time_scale$interval=="Ludlow"]]);
	fuzzy_ranges$region[!fuzzy_ranges$region %in% devonian_regions] <- "other";
	regional_richness <- hist(match(fuzzy_ranges$region[fuzzy_ranges$region!="other"],devonian_regions),breaks=0:length(devonian_regions),plot=F)$counts;
	names(regional_richness) <- devonian_regions;
	regional_richness <- sort(regional_richness,decreasing = TRUE);
	if (length(halloween_colors)<length(devonian_regions))	{
		group_colors <- c(halloween_colors,"magenta","cyan","blue","springgreen","lightslateblue")[1:length(devonian_regions)];
		group_colors <- c(group_colors,"grey");
		names(group_colors) <- c(names(regional_richness),"other");
		otu_information$otu_colors <- otu_colors <- group_colors[match(fuzzy_ranges$region,names(group_colors))];
		}
  }
taxon_group <- unique(otu_information[,plot_col]);
#otu_information[otu_information[,plot_col]=="" & otu_information$taxon_type=="ingroup",plot_col] <- "Trochocystitidae";
outgroup_ht <- sort(unique(otu_information[tolower(otu_information$taxon_type)=="outgroup",plot_col]));
ingroup_ht <- taxon_group[!taxon_group %in% outgroup_ht];
group_colors <- c();
nt <- 0;
noutgr <- length(outgroup_ht);
if (split_outgroups)	{
	while (nt < noutgr)	{
		nt <- nt+1;
		group_colors <- c(group_colors,paste("gray",round(nt/(noutgr+1)*100),sep=""));
		}
	} else if (noutgr>0)	{
	nt <- nt+1;
	group_colors <- rep("gray75",noutgr);
	}
ningr <- length(ingroup_ht);
#group_colors <- c(group_colors,rainbow(round((ningr/0.85),0))[1:ningr]);
if (length(halloween_colors)>length(ingroup_ht))	{
	group_colors <- c(group_colors,halloween_colors[1:length(ingroup_ht)]);
	} else	{
	group_colors <- c(group_colors,rainbow(length(ingroup_ht)+1)[1:length(ingroup_ht)])
	}
#group_colors <- c("gray75","red","orange","yellow")
ingroup_ht[ingroup_ht==""] <- ingroup;
names(group_colors) <- c(outgroup_ht,ingroup_ht);
group_colors <- group_colors[match(taxon_group,names(group_colors))]
otu_colors <- group_colors[match(otu_information[,plot_col],names(group_colors))];
names(otu_colors) <- otu_information$taxon;
otu_colors[is.na(otu_colors)] <- "black";

# rewritten 2020-12-16 to use hpd_age instead of branch_durations, which don't seem to add up.
root_ma <- min(hpd_age$md);
branch_ranges <- data.frame(start=as.numeric(vector_tree),
														finish=as.numeric(vector_tree));
rownames(branch_ranges) <- names(branch_durations)[!is.na(names(branch_durations))];
branch_ranges$finish <- hpd_age$md;
#cbind(branch_ranges$finish[1:notu],-fuzzy_ranges$fa_lb,-fuzzy_ranges$fa_ub)
#max(-fuzzy_ranges$fa_lb-branch_ranges$finish[1:notu])
#otu_names
branch_ranges[notu+1,] <- root_ma;
branch_ranges$start[1:notu] <- hpd_age$md[vector_tree[1:notu]]; # get divergence time of node
branch_ranges$start[(notu+2):(notu+nNodes)] <- hpd_age$md[vector_tree[(notu+2):(notu+nNodes)]];
problem_taxa <- (1:notu)[branch_ranges$start[1:notu]>branch_ranges$finish[1:notu]];

# double check branch_ranges
#branch_ranges$finish[(1:notu)[!(1:notu) %in% problem_taxa]];
#-abs(fuzzy_ranges$fa_lb[(1:notu)[!(1:notu) %in% problem_taxa]])-branch_ranges$finish[(1:notu)[!(1:notu) %in% problem_taxa]];

# get origin time for branches: if these are after last certain appearance of ancestor, then call it anagenetic
branch_origins <- branch_ranges$start;
branch_origins[1:notu] <- branch_ranges$finish[1:notu];
ancestral_spc <- (1:notu)[sampled_ancestors==1];
names(ancestral_spc) <- names(sampled_ancestors)[sampled_ancestors==1];
anagenetic_ancestors <- rep(0,notu);
names(anagenetic_ancestors) <- names(sampled_ancestors);
min_ranges$la_mr <- -abs(min_ranges$la_mr); min_ranges$fa_mr <- -abs(min_ranges$fa_mr);
anagenetic_sets <- data.frame(ancestor=as.character(),descendant=as.character());
anagenetic_sets_no <- data.frame(ancestor=as.numeric(),descendant=as.numeric());
# heavily rewritten 2023-10-09
if (length(ancestral_spc)>0)	{
#	pbdb_finds <- pbdb_data_list$pbdb_finds;
#	pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
	venn_tree <- transform_vector_tree_to_venn_tree(vector_tree);
	mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
#	species_origins <- branch_ranges$finish[1:notu];
#	cinctan_finds <- read.csv(paste(this_directory,"Cincta_Finds.csv",sep=""),header=TRUE,stringsAsFactors = FALSE);
	for (an in 1:length(ancestral_spc))	{
		anc <- ancestral_spc[an];
		daughter_node <- max(which(venn_tree==anc,arr.ind = TRUE)[,1]);
		its_node <- which(mat_tree==anc,arr.ind = TRUE)[,1]
		desc <- mat_tree[its_node,mat_tree[its_node,]!=anc]; # get descendants
	#	hpd_age[anc,]
	#	hpd_age[desc,]
	#	fuzzy_ranges[anc,]
	#	fuzzy_ranges[desc,]
		# prior routine using venn_tree
		progeny <- venn_tree[daughter_node,venn_tree[daughter_node,]>0];
		progeny <- progeny[progeny!=anc];
		# there can be only one anagenetic descendant: the last one. 
		if (min_ranges$la_mr[anc] < max(branch_origins[desc]) & min_ranges$la_mr[anc]!=0)  {
      anagenetic_ancestors[anc] <- 1;
  		} else if (min_ranges$la_mr[anc]==0 & -abs(fuzzy_ranges$la_lb[anc]) < max(branch_origins[desc])) {
      anagenetic_ancestors[anc] <- 1;
  		}
		if (anagenetic_ancestors[anc]==1) {
      # get the youngest descendant if there are 2+
		  desc <- desc[branch_origins[desc]==max(branch_origins[desc])];
		  dummy_set <- data.frame(ancestor=rownames(branch_ranges)[anc],
		                          descendant=rownames(branch_ranges)[desc]);
      dummy_set_no <- data.frame(ancestor=anc,descendant=desc);
      anagenetic_sets <- rbind(anagenetic_sets,dummy_set);
      anagenetic_sets_no <- rbind(anagenetic_sets_no,dummy_set_no);
      }
#		anc_find_data <- subset(pbdb_finds,pbdb_finds$accepted_name==otu_names[anc])
#		if(sum(abs(fuzzy_ranges$la_lb[anc])<=abs(fuzzy_ranges$fa_ub[progeny]))==0)	{
#			find_data <- subset(pbdb_finds,pbdb_finds$accepted_name %in% otu_names[c(anc,progeny)]);
#			accepted_names <- unique(pbdb_taxonomy$accepted_name[pbdb_taxonomy$taxon_name %in% otu_names[c(anc,progeny)]]);
#			accepted_names <- accepted_names[!is.na(accepted_names)];
#			poss_names <- pbdb_taxonomy$taxon_name[pbdb_taxonomy$accepted_name %in% accepted_names];
#			poss_taxon_nos <- pbdb_taxonomy$taxon_no[pbdb_taxonomy$accepted_name %in% accepted_names];
#			find_data <- pbdb_finds[pbdb_finds$accepted_no %in% poss_taxon_nos,]
#			find_data$accepted_name <- pbdb_taxonomy$accepted_name[match(find_data$accepted_no,pbdb_taxonomy$taxon_no)];
#			pbdb_finds[pbdb_finds$identified_name=="Crassatella tumidula",]
#			coccur_mat <- accersi_cooccurence_matrix(find_data);
#			if (sum(coccur_mat[match(otu_names[anc],rownames(coccur_mat)),])==1)
#			if (coccur_mat[1,2]==0)
#				anagenetic_ancestors[anc] <- 1;
#			}
		}
	}

#nNodes <- nrow(mat_tree)
phylo_axis <- get_phylo_axis_from_newick_string_w_anagenesis(newick_string = tree_output$newick_modified,sampled_ancestors=sampled_ancestors,anagenetic_ancestors=anagenetic_ancestors);
names(phylo_axis) <- c((otu_names),paste("node_",1:nNodes,sep=""));
if (notu>9)	names(phylo_axis)[notu+1:9] <- paste("node_0",1:9,sep="");

### use this to identify which sequences get lumped because of ancestry ####
# get anagenetic series of species
anagenetic_anc_no <- (1:notu)[anagenetic_ancestors==1];
names(anagenetic_anc_no) <- names(anagenetic_ancestors)[anagenetic_anc_no];
axis_ranks <- rank(phylo_axis[1:notu]);
axis_ranks <- axis_ranks-min(axis_ranks)+1;
anagenetic_overlaps <- hist(axis_ranks[1:notu],breaks=(min(axis_ranks[1:notu])-1):ceiling(max(axis_ranks[1:notu])),plot=FALSE)$counts;
#anagenetic_overlaps <- hist(phylo_axis[1:notu],breaks=((min(phylo_axis[1:notu])-1):max(phylo_axis[1:notu])),plot=FALSE)$counts;
anagenetic_segments <- names(anagenetic_ancestors[anagenetic_ancestors==1]);
#anagenetic_segments <- otu_names[anagenetic_ancestors==1];
anagenetic_segments_nex <- gsub(" ","_",anagenetic_segments);
#anagenetic_segments_nex <- otu_names_nex[anagenetic_ancestors==1];
names(phylo_axis) <- all_names <- rownames(branch_ranges);
names(anagenetic_overlaps) <- sort(unique(phylo_axis[1:notu]));
anagenetic_overlaps <- anagenetic_overlaps[anagenetic_overlaps>1];
anagenetic_fuzzies <- data.frame(ancestor=as.character(),descendant=as.character(),
								 fa_lb=as.numeric(),fa_ub=as.numeric(),la_lb=as.numeric(),la_ub=as.numeric(),
								 anc_no=as.numeric(),dsc_no=as.numeric(),
								 phylo_axis=as.numeric(),taxon=as.character(),stringsAsFactors = FALSE);
an <- 0;
while (an<length(anagenetic_overlaps))	{
	an <- an+1;
	dummy <- data.frame(ancestor=as.character("a"),descendant=as.character("z"),
						fa_lb=as.numeric(0),fa_ub=as.numeric(0),la_lb=as.numeric(0),la_ub=as.numeric(0),
						anc_no=as.numeric(0),dsc_no=as.numeric(0),phylo_axis=as.numeric(0),taxon=as.character("tx"),stringsAsFactors = FALSE);
	anagenetic_fuzzies <- rbind(anagenetic_fuzzies,dummy);
	anagenetic_fuzzies$ancestor[an] <- anagenetic_sets$ancestor[an];
	anagenetic_fuzzies$descendant[an] <- anagenetic_sets$descendant[an];
	morphospc <- (1:notu)[phylo_axis[1:notu]==as.numeric(names(anagenetic_overlaps)[an])];
#	if (length(unique(fuzzy_ranges$fa_lb[morphospc]))==1 && length(unique(fuzzy_ranges$fa_ub[morphospc])) && length(unique(fuzzy_ranges$la_lb[morphospc])))	{
	anc_no <- anagenetic_sets_no$ancestor[an];
	dsc_no <- anagenetic_sets_no$descendant[an];
#	anc_no <- morphospc[morphospc %in% anagenetic_anc_no];
#	dsc_no <- morphospc[!morphospc %in% anagenetic_anc_no];
	#	all_names[morphospc][!all_names[morphospc] %in% anagenetic_segments_nex]
#	anagenetic_fuzzies$ancestor[an] <- all_names[anc_no];
#	anagenetic_fuzzies$descendant[an] <- all_names[dsc_no];
	anagenetic_fuzzies$fa_lb[an] <- fuzzy_ranges$fa_lb[anc_no];
	anagenetic_fuzzies$fa_ub[an] <- max(fuzzy_ranges$fa_ub[morphospc]);
	anagenetic_fuzzies$la_lb[an] <- min(fuzzy_ranges$la_lb[morphospc]);
	anagenetic_fuzzies$la_ub[an] <- fuzzy_ranges$la_ub[dsc_no];
	anagenetic_fuzzies$anc_no[an] <- anc_no;
	anagenetic_fuzzies$dsc_no[an] <- dsc_no;
	anagenetic_fuzzies$taxon[an] <- paste(otu_names[anc_no],"->",otu_names[dsc_no],sep="");
	anagenetic_fuzzies$phylo_axis[an] <- phylo_axis[anc_no];
#		}
	}

#### plot time scale ####
st <- sum(max(abs(branch_ranges$start)) <= abs(finest_chronostrat$ma_lb));	# starts!
en <- 1+ sum(min(abs(fuzzy_ranges$la_ub)) < abs(finest_chronostrat$ma_ub));	# endsd!

standard_time_scale <- finest_chronostrat[st:en,];
standard_time_scale$ma_lb <- -abs(standard_time_scale$ma_lb);
standard_time_scale$ma_ub <- -abs(standard_time_scale$ma_ub);
standard_time_scale <- standard_time_scale[standard_time_scale$interval_sr=="",]

# Add names/symbols for plotting
strat_names <- standard_time_scale$st;
# Get the colors for time units
strat_colors <- standard_time_scale$color;
# Set the oldest and youngest intervals by names on whatever chronostratigraphic scale you used in the analysis
oldest_interval <- standard_time_scale$interval[1];
youngest_interval <- standard_time_scale$interval[nrow(standard_time_scale)];
# get the time scale that will be plotted
time_scale_to_plot <- unique(c(standard_time_scale$ma_lb,standard_time_scale$ma_ub));
# here, I want to adjust the earliest date: I don't want the whole Cambrian
yearbreaks <- as.numeric(set_axis_breaks_new(max(abs(branch_ranges$start))-min(fuzzy_ranges$la_ub)))
time_scale_to_plot[1] <- -ceiling(max(abs(branch_ranges$start))/max(yearbreaks))*max(yearbreaks);
#time_scale_to_plot[1] <- -abs(standard_time_scale$ma_lb[sum(-abs(standard_time_scale$ma_lb) < (-ceiling(max(abs(branch_ranges$start)))))]);

# Now, set the onset and end of the x-axis
onset <- min(time_scale_to_plot);
end <- max(time_scale_to_plot);

use_strat_labels <- TRUE;				# if TRUE, then strat_names will be plotted on X-axis inside boxes
alt_back <- FALSE;								# if TRUE, then the background will alternat shades between major intervals
plot_title <- "";							# Name of the plot; enter "" for nothing
ordinate <- "";								# Label of Y-axis
hues <- TRUE;									# If TRUE, then IGN stratigraphic colors will be used
colored <- "base";							# Where IGN stratigraphic colors should go
end_extra <- abs(end-onset)/4;
ysize <-  min(6,4.285714285*notu/20);
#ysize <- 5;
xsize <- 1.25*5; 
mxy <- max(phylo_axis)
mny <- min(phylo_axis);									# set maximum y value
ysize <- 5;
print_names <- TRUE;
print_pix <- FALSE;
strat_label_size <- 1;					# size of the labels for chronostratigraphic units
myr_size <- 1;
Phanerozoic_Timescale_Plot_Flexible(onset,end+end_extra,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=2*strat_label_size,font=franky,xlab_size=1.5,myr_size=myr_size);
#### Plot Tree ####
sp <- 0;
#otu_colors[otu_colors=="white"] <- "gray75";
while (sp < notu) {
#for (sp in 1:notu)	{
	sp <- sp+1;
	if (!sp %in% anagenetic_fuzzies$anc_no & !sp %in% anagenetic_fuzzies$dsc_no)	{
		# upper limit exceeds lower limit, but no definite range
		if (fuzzy_ranges$fa_lb[sp]==fuzzy_ranges$la_lb[sp] && fuzzy_ranges$fa_ub[sp]>fuzzy_ranges$la_ub[sp])	{
			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=50),border=makeTransparent(otu_colors[sp],alpha=50),lwd=0);
#			print(sp);
			} else if (fuzzy_ranges$fa_ub[sp]==fuzzy_ranges$la_ub[sp] && fuzzy_ranges$fa_lb[sp]>fuzzy_ranges$la_lb[sp])	{
		# lower limit exceeds exceeds limit, but no definite range
			rect(-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=50),border=makeTransparent(otu_colors[sp],alpha=50),lwd=0);
#			print(sp);
			} else	{
			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			}
		if (abs(fuzzy_ranges$fa_ub[sp])>abs(fuzzy_ranges$la_lb[sp]))
			rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+1/3,col=otu_colors[sp],border=otu_colors[sp],lwd=0);
#		if (print_names)	text(-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp],fuzzy_ranges$taxon[sp],family=franky,pos=4,cex=0.5);
#		if ((fuzzy_ranges$fa_lb[sp]==fuzzy_ranges$la_lb[sp]) || (fuzzy_ranges$fa_ub[sp]==fuzzy_ranges$la_ub[sp]))	{
#			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100));
#			} else	{
#			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100));
#			rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+1/3,col=otu_colors[sp],border=otu_colors[sp]);
#			}
		} else if (sp %in% anagenetic_fuzzies$dsc_no)	{
		fuzzy_dude <- match(sp,anagenetic_fuzzies$dsc_no);
		rect(-abs(anagenetic_fuzzies$fa_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]-1/3,-abs(anagenetic_fuzzies$la_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
		if (abs(anagenetic_fuzzies$fa_ub[fuzzy_dude])>abs(anagenetic_fuzzies$la_lb[fuzzy_dude]))
			rect(-abs(anagenetic_fuzzies$fa_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]-1/3,-abs(anagenetic_fuzzies$la_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+1/3,col=otu_colors[sp],lwd=0);
#		if ()
		}
	}
observed_nodes_all <- rep(0,nNodes);
an <- 0;
while (an < length(ancestral_spc))	{
	an <- an+1;
	observed_nodes_all[which(mat_tree==ancestral_spc[an],arr.ind = TRUE)[1]] <- ancestral_spc[an];
	}
for (nn in 1:length(phylo_axis))	{
	if (nn<=notu && !nn %in% problem_taxa)	{
		segments(branch_ranges$start[nn],phylo_axis[nn],branch_ranges$finish[nn],phylo_axis[nn],lwd=1,lty=2)
		} else if (nn>(notu+1))	{
		nd <- nn-notu;
		if (observed_nodes_all[nd]==0)	{
			segments(branch_ranges$start[nn],phylo_axis[nn],branch_ranges$finish[nn],phylo_axis[nn],lwd=5*clade_posteriors[nd]);
			} else	{
			segments(branch_ranges$start[nn],phylo_axis[nn],branch_ranges$finish[nn],phylo_axis[nn],lwd=5*clade_posteriors[nd],lty=3);
			}
		} # end case for clade
	}
for (nd in 1:nNodes)	{
	htu <- nd+notu;
	f1 <- mat_tree[nd,mat_tree[nd,]>0];
	segments(branch_ranges$finish[htu],min(phylo_axis[f1]),branch_ranges$finish[htu],max(phylo_axis[f1]));
	}
# plot inferred origin times of otus
for (sp in 1:notu)	points(-abs(branch_ranges$finish[sp]),phylo_axis[sp],pch=8,cex=10/(notu/2));
for (sp in 1:notu)	points(-abs(branch_ranges$finish[sp]),phylo_axis[sp],pch=8,cex=10/(notu/2),lwd=0.25,col=otu_colors[sp]);
top <- mxy+1;
otu_font_size <- max(0.5,ysize*25/(6*notu));
if (print_names)	{
#	otu_names <- gsub("Crassatellites","Crassatella",otu_names);
	sp <- 0;
#	for (sp in 1:notu)	{
	while (sp < notu)	{
		sp <- sp+1;
		if (!sp %in% anagenetic_fuzzies$anc_no && !sp %in% anagenetic_fuzzies$dsc_no)	{
			xx <- max(branch_ranges$finish[sp],-abs(fuzzy_ranges$la_ub[sp]));
			text(xx,phylo_axis[sp],otu_names[sp],pos=4,cex=otu_font_size,family=franky);
			} else	if (sp %in% anagenetic_fuzzies$dsc_no)	{
			fuzzy_dude <- match(sp,anagenetic_fuzzies$dsc_no);
			if (is.na(fuzzy_dude))
			  fuzzy_dude <- match(sp,anagenetic_fuzzies$anc_no);
#		lineage_name <- paste(anagenetic_fuzzies$ancestor[fuzzy_dude],"->",anagenetic_fuzzies$descendant[fuzzy_dude],sep="")
			xx <- max(branch_ranges$finish[sp],-abs(fuzzy_ranges$la_ub[sp]));
			text(xx,phylo_axis[sp],anagenetic_fuzzies$taxon[fuzzy_dude],pos=4,cex=otu_font_size,family=franky);
			}
		}
	}
# ####
if (print_pix)	{
	pictured <- c("Protocinctus_mansillaensis","Elliptocinctus_barrandei","Trochocystites_bohemicus");
	x_add <- abs(onset-end)/100;
	for (pc in 1:length(pictured))	{
		pic_file <- paste(this_directory,pictured[pc],"_horizontal.png",sep="");
		png_info <- readPNG(pic_file);
		txn <- match(pictured[pc],otu_names_nex);
		add_png(png_info, x = -abs(fuzzy_ranges$la_ub[txn])+x_add,y=phylo_axis[txn],height = 1,x_cent=FALSE)
		}
	}
key_top <- TRUE;
onset_right <- (end-onset)/10;
onset_right2 <- 19*(end-onset)/200;
group_font_size <- min(2*otu_font_size,1.5);
higher_taxon_key_colors <- group_colors[match(c(outgroup_ht,ingroup_ht),names(group_colors))];
if (length(outgroup_ht)>1 && !split_outgroups)	{
	composite_outgroup <- paste(outgroup_ht,collapse=" + ");
	} else	{
	composite_outgroup <- outgroup_ht;
	}
names(higher_taxon_key_colors)[names(higher_taxon_key_colors) %in% outgroup_ht] <- composite_outgroup;
xx <- unique(names(higher_taxon_key_colors));
higher_taxon_key_colors <- unique(higher_taxon_key_colors);
names(higher_taxon_key_colors) <- xx;
dt <- notu*0.025;
if (!key_top)	{
	top <- 1.25*dt*length(higher_taxon_key_colors);
	} else	{
	top <- notu;
	}
for (ii in 1:length(higher_taxon_key_colors))	{
#	rect(onset,top,onset+onset_right,top-1);
#	makeTransparent(higher_taxon_key_colors[ii],alpha=100)
	rect(onset,top,onset+onset_right/2,top-1,col=makeTransparent(higher_taxon_key_colors[ii],alpha=100),lwd=0);
	rect((onset+onset_right/2),top,onset+onset_right,top-1,col=higher_taxon_key_colors[ii],lwd=0);
	rect(onset,top,onset+onset_right,top-1);
	text(onset+onset_right2,top-0.5,names(higher_taxon_key_colors)[ii],pos=4,family=franky,cex=group_font_size);
	top <- top-1.25*dt;
	} 

if (!key_top)	{
	top <- 2*1.25*dt*length(higher_taxon_key_colors);
	} else	{
	top <- notu-2.5*1.25*dt*length(higher_taxon_key_colors);
	}

post_p <- 0.875;
for (i in 1:4)	{
	segments(onset,top,onset+onset_right,lwd=5*post_p);
	text(onset+onset_right2,top,paste("posterior prob=",round(post_p,3)),pos=4,family=franky,cex=group_font_size);
	top <- top-1.25*dt;
	post_p <- post_p-0.25;
	}

names(vector_tree) <- names(phylo_axis);
cinctan_family_nodes <- vector(length=length(higher_taxa));
names(cinctan_family_nodes) <- higher_taxa;
names(cinctan_family_nodes)[1] <- "Cincta"
cinctan_family_nodes[1] <- 25;
cinctan_family_nodes[2] <- 38;
cinctan_family_nodes[3] <- 30;
cinctan_family_nodes[4] <- 27;

for (i in 1:4)	{
	cfn <- cinctan_family_nodes[i];
#	rect(hpd_age$lb[cfn],phylo_axis[cfn]-0.125,hpd_age$ub[cfn],phylo_axis[cfn]+0.125,col="white",lwd=0);
#	if (i==1)	{
#		temp_col <- makeTransparent("magenta",alpha=50);
#		rect(hpd_age$lb[cfn],phylo_axis[cfn]-0.125,hpd_age$ub[cfn],phylo_axis[cfn]+0.125,col=temp_col,lwd=0);
#		} else	{
		temp_col <- makeTransparent(higher_taxon_colors[i],alpha=100);
		rect(hpd_age$lb[cfn],phylo_axis[cfn]-0.125,hpd_age$ub[cfn],phylo_axis[cfn]+0.125,col=temp_col,lwd=0);
#		}
	}
#for (pc in 1:length(pictured))	{
#	jpg_file <- paste(this_directory,pictured[pc],"_horizontal.jpg",sep="");
#	txn <- match(pictured[pc],otu_names_nex);
#	add_jpeg(p,x=-abs(fuzzy_ranges$la_ub[txn])+x_add,y=phylo_axis[txn]-0.5,add=TRUE)
#	}

