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

options(warn=0);
load(paste(data_for_R_folder,"Gradstein_2020_Augmented.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
load(paste(data_for_R_folder,"Paleobiology_Database.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
load(paste(data_for_R_folder,"PaleoDB_Edits.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
load(paste(data_for_R_folder,"Rock_Unit_Database.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
options(warn=1);

pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
pbdb_taxonomy <- rbind(pbdb_taxonomy,paleodb_fixes$paleodb_taxonomy_edits[!paleodb_fixes$paleodb_taxonomy_edits$taxon_no %in% pbdb_taxonomy$taxon_no,]);
pbdb_finds <- pbdb_data_list$pbdb_finds;
pbdb_sites <- pbdb_data_list$pbdb_sites_refined;
#pbdb_sites$formation[pbdb_sites$collection_no==210883] <- "Ihandar";
#pbdb_sites$rock_no[pbdb_sites$collection_no==210883] <- pbdb_sites$rock_no_sr[pbdb_sites$collection_no==210883] <- 8329;
rock_database <- rock_unit_data$rock_unit_database;

gap <- INAP <- -22;
UNKNOWN <- missing <- -11;
polymorphs <- TRUE	# if false, then these are converted to unknowns
split_outgroups <- FALSE;

# set parameters & files for printing tree ####
variable_to_plot <- "taxonomy"; # use "taxonomy","geography","ecology"
variable_name <- "higher_taxon";
rank_to_plot <- "subfamily";
organize_by_parent_taxon <- TRUE;
strat_rank <- "Stage";		# I recommend either "Stage" or "Substage"
stratigraphic_scale <- "International";		# I recommend either "International" or "Stage Slice"
mctree_file <- paste(tree_drawing_folder,"max_clade_cred_tree_FBDMM.tre",sep="");
#stratigraphic_range_file <- paste(tree_drawing_folder,"Acanthoceratidae_Range_Data.xlsx",sep="");
elevate_subspecies <- FALSE;

#### Get Tree Information for Plotting ####
mctree <- scan(file=mctree_file,what=character(),sep="\n");
#mctree <- scan(file=file.choose(),what=character(),sep="\n");
mctree <- gsub("\t","",mctree);
begin_taxa <- match("begin taxa;",tolower(mctree));
taxlabels <- match("taxlabels",tolower(mctree));
notu <- as.numeric(strsplit(gsub(";","",mctree[begin_taxa+1]),split = "=")[[1]][2]);
otu_names_nex <- mctree[taxlabels+(1:notu)];
otu_names <- gsub("_"," ",otu_names_nex);
newick_string_full <- mctree[length(mctree)-1];
newick_string_full <- paste(strsplit(newick_string_full,"")[[1]][match("(",strsplit(newick_string_full,"")[[1]]):length(strsplit(newick_string_full,"")[[1]])],collapse="");
if (gsub(otu_names_nex[1],"",newick_string_full)!=newick_string_full)	{}
mctree_info <- read_newick_string_mcmc_numbered(newick_string_full,otu_names);
branch_durations <- mctree_info$branch_durations;
clade_posteriors <- mctree_info$clade_posteriors;
sampled_ancestors <- mctree_info$ancestral;
vector_tree <- mctree_info$vector_tree;
nNodes <- length(vector_tree)-notu;
htu_heights <- accersi_patristic_distance_from_base(vector_tree);
br_lngth_history <- array(0,dim=c(notu+nNodes,1+max(htu_heights)));
rownames(br_lngth_history) <- names(vector_tree);
hpd_age <- mctree_info$hpd_age;
for (sp in 1:(notu+nNodes))	{
	anc_node <- vector_tree[sp];
	pd <- 1;
	br_lngth_history[sp,pd] <- branch_durations[sp];
#	pd <- pd+1;
#	br_lngth_history[sp][pd] <- branch_durations[anc_node];
	while (anc_node!=-1)	{
		pd <- pd+1;
		br_lngth_history[sp,pd] <- branch_durations[anc_node];
		anc_node <- vector_tree[anc_node];
		}
	}
br_lngth_history <- br_lngth_history[,colSums(br_lngth_history)>0];
estimated_heights <- rowSums(br_lngth_history);
estimated_depths <- -(estimated_heights-max(estimated_heights));
#cbind(hpd_age,estimated_heights,estimated_depths);
# get taxonomic information ####
otu_taxonomy <- pbdb_taxonomy[match(otu_names,pbdb_taxonomy$taxon_name)[!is.na(match(otu_names,pbdb_taxonomy$taxon_name))],];
if (nrow(otu_taxonomy)==0)	{
	# Added to deal with FBDMM trees 2024-05-22
#	otu_names <- gsub("_"," ",otu_names_nex);
	genus_names <- sapply(otu_names,divido_genus_names_from_species_names);
	taxon_ranks <- pbdb_taxonomy$accepted_rank[match(genus_names,pbdb_taxonomy$taxon_name)];
	if (sum(taxon_ranks %in% "species")>0)	{
		otu_names_orig <- otu_names;
		addenda <- array("",dim=c(notu,1));
		for (sp in 1:notu)	{
			if (!is.na(taxon_ranks[sp]) & taxon_ranks[sp]=="species")	{
				extra <- gsub(paste(genus_names[sp],""),"",otu_names[sp]);
				if (length(strsplit(extra," ")[[1]])>ncol(addenda))	{
					a <- length(strsplit(extra," ")[[1]])-ncol(addenda);
					dummy <- array("",dim=c(notu,a));
					addenda <- cbind(addenda,dummy);
					}
				addenda[sp,] <- strsplit(extra," ")[[1]];
				otu_names[sp] <- genus_names[sp];
				}
			}
		if (sum(!taxon_ranks %in% c("genus","subgenus"))>0)	{
			for (sp in 1:notu)	{
				if (is.na(taxon_ranks[sp]) || taxon_ranks[sp] %in% c("genus","subgenus"))	{
					breakdown <- strsplit(otu_names[sp]," ")[[1]];
					b <- length(breakdown);
					addenda[sp,] <- breakdown[(b-1):b];
					otu_names[sp] <- paste(breakdown[1:(b-2)],collapse=" ");
					}
				}
			}
		}
	otu_taxonomy <- pbdb_taxonomy[match(otu_names,pbdb_taxonomy$taxon_name)[!is.na(match(otu_names,pbdb_taxonomy$taxon_name))],];
	}

missing_species <- otu_names[is.na(match(otu_names,pbdb_taxonomy$taxon_name))];
ms <- 0;
while(ms < length(missing_species))	{
	ms <- ms+1;
	genus <- divido_genus_names_from_species_names(missing_species[ms]);
	genus_members <- pbdb_taxonomy[pbdb_taxonomy$genus %in% genus,];
	if (nrow(genus_members)==0)	{
		new_genus <- pbdb_taxonomy$accepted_name[pbdb_taxonomy$taxon_name==genus];
		if (length(new_genus)==0)	{
			genus_members <- pbdb_taxonomy[gsub(genus,"",pbdb_taxonomy$taxon_name)!=pbdb_taxonomy$taxon_name,];
#			genus_members <- genus_members[genus_members$taxon_rank %in% c("species","subspecies"),];
			} else	{
#			pbdb_taxonomy[pbdb_taxonomy$taxon_name==new_genus,]
			genus_members <- pbdb_taxonomy[pbdb_taxonomy$parent_name %in% c(new_genus,genus),];
			}
#		genus_members <- pbdb_taxonomy[gsub(genus,"",pbdb_taxonomy$genus)!=pbdb_taxonomy$genus,];
		}
	genus_members <- accersi_faux_genus_species_combos_for_species_in_subgenera(paleodb_taxonomy=genus_members);
	if (!is.na(match(missing_species[ms],genus_members$taxon_name)))	{
		otu_taxonomy <- rbind(otu_taxonomy,genus_members[match(missing_species[ms],genus_members$taxon_name),]);
		pbdb_taxonomy <- rbind(pbdb_taxonomy,genus_members[match(missing_species[ms],genus_members$taxon_name),]);
		}
	}
missing_species <- otu_names[is.na(match(otu_names,pbdb_taxonomy$taxon_name))];
higher_taxon_no <- accersi_common_higher_taxon_no(pbdb_taxonomic_data=otu_taxonomy,pbdb_taxonomy);
higher_taxon <- pbdb_taxonomy$taxon_name[pbdb_taxonomy$taxon_no==higher_taxon_no];
higher_taxon_rank <- pbdb_taxonomy$taxon_rank[pbdb_taxonomy$taxon_no==higher_taxon_no];
higher_taxon_rank_no <- paste(higher_taxon_rank,"_no",sep="");
if (!higher_taxon_rank %in% colnames(pbdb_taxonomy))	{
	highest_present <- standard_pbdb_taxon_ranks[match(standard_pbdb_taxon_ranks,taxonomic_rank)>match(higher_taxon_rank,taxonomic_rank)][length(standard_pbdb_taxon_ranks[match(standard_pbdb_taxon_ranks,taxonomic_rank)>match(higher_taxon_rank,taxonomic_rank)])];
	highest_present_no <- paste(highest_present,"_no",sep="");
	ingroup_taxonomy <- pbdb_taxonomy[pbdb_taxonomy[,highest_present_no] %in% otu_taxonomy[,highest_present_no][1],];
#	otu_names[!otu_names %in% group_taxonomy$taxon_name]
	ingroup_taxonomy <- add_higher_taxon_to_pbdb_taxonomy(taxon_rank=higher_taxon_rank,paleodb_taxonomy=ingroup_taxonomy);
#	otu_names[!otu_names %in% group_taxonomy$taxon_name]
	} else	{
	ingroup_taxonomy <- pbdb_taxonomy[pbdb_taxonomy[,higher_taxon_rank_no] %in% otu_taxonomy[,higher_taxon_rank_no][1],];
	}
ingroup_taxonomy <- ingroup_taxonomy[ingroup_taxonomy[,higher_taxon_rank] %in% higher_taxon,];
ingroup_taxonomy <- pretend_subgenera_are_genera(pbdb_taxonomic_data = ingroup_taxonomy);

if (!rank_to_plot %in% colnames(ingroup_taxonomy))
	ingroup_taxonomy <- add_higher_taxon_to_pbdb_taxonomy(taxon_rank=rank_to_plot,paleodb_taxonomy=ingroup_taxonomy);

fuzzy_ranges <- data.frame(taxon=otu_names,fa_lb=as.numeric(rep(0,notu)),fa_ub=as.numeric(rep(0,notu)),la_lb=as.numeric(rep(0,notu)),la_ub=as.numeric(rep(0,notu)));
for (sp in 1:notu)	{
	otu_nos <- c(ingroup_taxonomy$taxon_no[match(fuzzy_ranges$taxon[sp],ingroup_taxonomy$taxon_name)],ingroup_taxonomy$accepted_no[match(fuzzy_ranges$taxon[sp],ingroup_taxonomy$taxon_name)]);
	otu_finds <- pbdb_finds[pbdb_finds$accepted_no %in% otu_nos,];
	otu_sites <- pbdb_sites[pbdb_sites$collection_no %in% otu_finds$collection_no,];
	fuzzy_ranges$fa_lb[sp] <- max(otu_sites$ma_lb);
	fuzzy_ranges$fa_ub[sp] <- max(otu_sites$ma_ub);
	fuzzy_ranges$la_lb[sp] <- min(otu_sites$ma_lb);
	fuzzy_ranges$la_ub[sp] <- min(otu_sites$ma_ub);
	}

fuzzy_ranges$taxon <- gsub("_"," ",fuzzy_ranges$taxon);
notu <- nrow(fuzzy_ranges);
#fuzzy_ranges[,c("region","realm")];
otu_information <- data.frame(taxon=fuzzy_ranges$taxon);
otu_information$genus <- ingroup_taxonomy$genus[match(otu_information$taxon,ingroup_taxonomy$taxon_name)];
otu_information$family <- ingroup_taxonomy$family[match(otu_information$taxon,ingroup_taxonomy$taxon_name)];
otu_information$order <- ingroup_taxonomy$order[match(otu_information$taxon,ingroup_taxonomy$taxon_name)];
otu_information$class <- ingroup_taxonomy$class[match(otu_information$taxon,ingroup_taxonomy$taxon_name)];
if (variable_to_plot=="taxonomy" & !rank_to_plot %in% colnames(otu_information))	{
	otu_information <- cbind(otu_information,ingroup_taxonomy[match(otu_information$taxon,ingroup_taxonomy$taxon_name),rank_to_plot]);
	colnames(otu_information)[ncol(otu_information)] <- rank_to_plot;
	}
otu_information$taxon_type <- rep("ingroup",notu);
otu_information$subfamily[otu_information$subfamily==""] <- gsub("idae","inae",otu_information$family[otu_information$subfamily==""]);
# use this for anagenesis.
min_ranges <- data.frame(taxon=fuzzy_ranges$taxon,fa_mr=rep(0,notu),la_mr=rep(0,notu));
for (sp in 1:notu)  {
	if (abs(fuzzy_ranges$fa_ub[sp])>=abs(fuzzy_ranges$la_lb[sp]))  {
		min_ranges$fa_mr[sp] <- fuzzy_ranges$fa_ub[sp];
		min_ranges$la_mr[sp] <- fuzzy_ranges$la_lb[sp];
		}
	}

finest_chronostrat <- gradstein_2020_emended$time_scale[gradstein_2020_emended$time_scale$scale==stratigraphic_scale,];
finest_chronostrat <- finest_chronostrat[finest_chronostrat$chronostratigraphic_rank==strat_rank & finest_chronostrat$interval_sr=="",];
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];
finest_chronostrat <- finest_chronostrat[finest_chronostrat$ma_lb>=min(fuzzy_ranges$la_ub),];

# remove any taxa from ranges that are not in the tree;
fuzzy_ranges <- fuzzy_ranges[fuzzy_ranges$taxon %in% otu_names,];
otu_information <- otu_information[otu_information$taxon %in% otu_names,];

faux_recent <- min(abs(fuzzy_ranges$fa_ub));
estimated_depths_ma <- -faux_recent-estimated_depths;
oldest_fa <- -(max(fuzzy_ranges$fa_ub)+max(fuzzy_ranges$fa_lb))/2;

rescaled_fa_bounds <- fuzzy_ranges;
rescaled_fa_bounds$la_lb <- rescaled_fa_bounds$la_ub <- NULL;
rescaled_fa_bounds$fa_lb <- rescaled_fa_bounds$fa_lb - faux_recent;
rescaled_fa_bounds$fa_ub <- rescaled_fa_bounds$fa_ub - faux_recent;
#sort(rescaled_fa_bounds$taxon)
basal_taxon <- fuzzy_ranges$taxon[match(max(fuzzy_ranges$fa_lb),fuzzy_ranges$fa_lb)];
basal_date_raw <- max(((mctree_info$hpd$lb+mctree_info$hpd$ub)/2)[1:notu]);	# base from Bayesian tree before rescaling 0 to actual upper bound
basal_taxon_no <- match(basal_date_raw,((mctree_info$hpd$lb+mctree_info$hpd$ub)/2)[1:notu]);

hpd_age <- -(mctree_info$hpd_age+faux_recent); # upper & lower bound of FA;
hpd_age$md <- (hpd_age$lb+hpd_age$ub)/2;

clade_onset <- max(abs(hpd_age));
#clade_onset <- max(abs(branch_ranges));
finest_chronostrat <- finest_chronostrat[finest_chronostrat$ma_ub<clade_onset,];
finest_chronostrat$st[finest_chronostrat$st=="Co1"] <- "L.Con";

#branch_durations <- branch_durations2;

mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree=mctree_info$vector_tree);
nclades <- nrow(mat_tree);
venn_tree <- transform_vector_tree_to_venn_tree(vector_tree);
nNodes <- nrow(mat_tree);
lowest_taxon <- min(which(mat_tree %in% (1:notu),arr.ind = TRUE)); # get the OTU that is closest to the base of the tree
outgp <- mat_tree[lowest_taxon,mat_tree[lowest_taxon,]<=notu];
#-max(abs(hpd$md[outgp]))
branch_dates <- data.frame(start=as.numeric(rep(0,nclades+notu)),
												 finish=as.numeric(rep(0,nclades+notu)));
branch_dates$finish[notu+1] <- branch_durations[notu+1];
rownames(branch_dates) <- names(branch_durations);
for (cl in 1:nclades)	{
	htu <- cl+notu;
	desc <- mat_tree[cl,][mat_tree[cl,]>0];
	branch_dates$start[desc] <- branch_dates$finish[htu];
#	branch_dates$finish[desc] <- branch_dates$finish[htu]+branch_durations2[desc];
	branch_dates$finish[desc] <- branch_dates$finish[htu]+branch_durations[desc];
	}

comp_sp <- match(rownames(hpd_age)[1],rownames(branch_dates));
branch_ranges <- branch_dates+(hpd_age$md[1]-branch_dates$finish[comp_sp]);
min(fuzzy_ranges$fa_lb-abs(branch_ranges$finish[1:notu]));
write.csv(cbind(abs(branch_ranges[1:notu,]),fuzzy_ranges),"WTF.csv",row.names = F);
branch_adj <- max(0,max(fuzzy_ranges$fa_ub-abs(branch_ranges$finish[1:notu])));
branch_ranges <- branch_ranges-branch_adj;
#rownames(branch_ranges) <- rownames(hpd);
#### Plot like a Lannister ####
if (variable_to_plot=="taxonomy")  {
	plot_col <- match(tolower(rank_to_plot),tolower(colnames(otu_information)));
	if (is.na(plot_col))	{
		otu_information <- cbind(otu_information,ingroup_taxonomy[match(otu_information$taxon,ingroup_taxonomy$taxon_name),tolower(rank_to_plot)]);
		colnames(otu_information)[ncol(otu_information)] <- rank_to_plot;
		colranks <- match(colnames(otu_information),taxonomic_rank);
		colranks[is.na(colranks)][1] <- -1;
		colranks[is.na(colranks)] <- 100;
		otu_information <- otu_information[,order(colranks)];
		if (tolower(rank_to_plot)=="subfamily")
			otu_information$subfamily[otu_information$subfamily==""] <- gsub("idae","inae",otu_information$family[otu_information$subfamily==""]);
		plot_col <- match(tolower(rank_to_plot),tolower(colnames(otu_information)));
		}
	orig_ranks <- colnames(otu_information)[colnames(otu_information) %in% taxonomic_rank];
	orig_ranks <- orig_ranks[!orig_ranks %in% rank_to_plot];
	parent_rank <- orig_ranks[match(orig_ranks,taxonomic_rank)>match(rank_to_plot,taxonomic_rank)][1];
	
	taxon_group <- unique(otu_information[,plot_col]);
	atomized_parent <- strsplit(higher_taxon,"")[[1]];
	lparent <- length(atomized_parent);
	parent_matches <- c();
	for (tg in 1:length(taxon_group))	{
		atomized_group <- strsplit(taxon_group[tg],"")[[1]];
		ldaughter <- length(atomized_group);
		parent_matches <- c(parent_matches,sum(atomized_group[1:min(ldaughter,lparent)]==atomized_parent[1:min(ldaughter,lparent)]));
		}
	names(parent_matches) <- taxon_group;
	eponym <- taxon_group[match(max(parent_matches),parent_matches)];
	non_eponyms <- sort(taxon_group[!taxon_group %in% eponym]);
	taxon_group <- c(eponym,non_eponyms);
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

# organize (say) subfamilies within families.
	if (organize_by_parent_taxon)	{
		hierarchy <- unique(otu_information[,c(match(parent_rank,colnames(otu_information)),plot_col+(0:1))]);
#		if (!parent_rank %in% colnames(hierarchy))
#			hierarchy <- cbind(hierarchy,otu_information[,parent_rank]);
		hierarchy$family[hierarchy$subfamily=="Drevermanniinae"] <- "Phillipsiidae (?)"
		hierarchy <- hierarchy[order(hierarchy[,2],hierarchy[,1]),];
		ingroup_ht <- hierarchy[,1];
		if (nrow(hierarchy)<=length(halloween_colors))	{
			hierarchy$color <- halloween_colors[1:nrow(hierarchy)];
			} else {
			hierarchy$color <- c("red","orange","gold2","yellow3","forestgreen","cyan","deepskyblue2","blue","purple1","gray")[1:nrow(hierarchy)];
			}
		group_colors <- hierarchy$color;
		} else {
		if (length(taxon_group)<=length(halloween_colors))	{
			group_colors <- halloween_colors[1:length(taxon_group)];
			} else {
			group_colors <- c("red","orange","gold2","yellow3","forestgreen","cyan","deepskyblue2","blue","purple1","gray")[1:nrow(hierarchy)];
			}
		}
	
#group_colors <- c(group_colors,rainbow(round((ningr/0.85),0))[1:ningr]);
#if (length(halloween_colors)>length(ingroup_ht))	{
#	group_colors <- c(group_colors,halloween_colors[1:length(ingroup_ht)]);
#	} else	{
#	group_colors <- c(group_colors,rainbow(length(ingroup_ht)+1)[1:length(ingroup_ht)])
#	}
#ingroup_ht[ingroup_ht==""] <- ingroup;
	names(group_colors) <- c(outgroup_ht,ingroup_ht);
	group_colors <- group_colors[match(taxon_group,names(group_colors))]
	otu_information$otu_colors <- otu_colors <- group_colors[match(otu_information[,plot_col],names(group_colors))];
	if (organize_by_parent_taxon) otu_information$otu_colors <- otu_colors <- hierarchy$color[match(otu_information[,plot_col],hierarchy$subfamily)]
	names(otu_colors) <- otu_information$taxon;
	otu_colors[is.na(otu_colors)] <- "black";
	} else if (variable_to_plot=="geography")  {
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

# Set up coordinates for branches & taxa ####
# rewritten 2020-12-16 to use hpd instead of branch_durations, which don't seem to add up.
# rerewritten 2024-03-08 because BEAST trees hpd is not working for me...
root_ma <- min(hpd_age$md);
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
	#	hpd[anc,]
	#	hpd[desc,]
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
phylo_axis <- get_phylo_axis_from_newick_string_w_anagenesis(newick_string = mctree_info$newick,sampled_ancestors=sampled_ancestors,anagenetic_ancestors=anagenetic_ancestors);
#phylo_axis <- phylo_axis-4;
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
as <- 0;
anagenetic_fuzzies <- data.frame(ancestor=as.character(),descendant=as.character(),
																 fa_lb=as.numeric(),fa_ub=as.numeric(),la_lb=as.numeric(),la_ub=as.numeric(),
																 anc_no=as.numeric(),dsc_no=as.numeric(),
																 phylo_axis=as.numeric(),taxon=as.character(),stringsAsFactors = FALSE);
while (as < nrow(anagenetic_sets))	{
	as <- as+1;
	annie <- anagenetic_segments[as];
	descy <- names(phylo_axis)[phylo_axis==phylo_axis[match(annie,names(phylo_axis))]];
	descy <- descy[!descy %in% annie];
	descy <- descy[descy %in% otu_names_nex];
	dummy <- data.frame(ancestor=as.character(annie),descendant=as.character("z"),
						fa_lb=as.numeric(0),fa_ub=as.numeric(0),la_lb=as.numeric(0),la_ub=as.numeric(0),
						anc_no=as.numeric(0),dsc_no=as.numeric(0),phylo_axis=as.numeric(0),taxon=as.character("tx"),stringsAsFactors = FALSE);
#	morphospc <- (1:notu)[phylo_axis[1:notu]==as.numeric(names(anagenetic_overlaps)[as])];
	morphospc <- (1:notu)[phylo_axis[1:notu]==phylo_axis[match(annie,names(phylo_axis[1:notu]))]];
	dummy$ancestor <- anagenetic_sets$ancestor[as];
	dummy$descendant <- anagenetic_sets$descendant[as];
	dummy$anc_no <- anagenetic_sets_no$ancestor[as];
	dummy$dsc_no <- anagenetic_sets_no$descendant[as];
	dummy$fa_lb <- fuzzy_ranges$fa_lb[dummy$anc_no];
	dummy$fa_ub <- max(fuzzy_ranges$fa_ub[morphospc]);
	dummy$la_lb <- min(fuzzy_ranges$la_lb[morphospc]);
	dummy$la_ub <- fuzzy_ranges$la_ub[dummy$dsc_no];
	dummy$taxon <- paste(otu_names[dummy$anc_no],"->",otu_names[dummy$dsc_no],sep="");
	dummy$phylo_axis <- phylo_axis[dummy$anc_no];
	if (dummy$dsc_no<=notu)
		anagenetic_fuzzies <- rbind(anagenetic_fuzzies,dummy);
	}

# old, not-quite working version
names(anagenetic_overlaps) <- sort(unique(phylo_axis[1:notu]));
anagenetic_overlaps <- anagenetic_overlaps[anagenetic_overlaps>1];
anagenetic_fuzzies <- data.frame(ancestor=as.character(),descendant=as.character(),
																 fa_lb=as.numeric(),fa_ub=as.numeric(),la_lb=as.numeric(),la_ub=as.numeric(),
																 anc_no=as.numeric(),dsc_no=as.numeric(),
																 phylo_axis=as.numeric(),taxon=as.character(),stringsAsFactors = FALSE);
an <- MAXNO;
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
finest_chronostrat <- gradstein_2020_emended$time_scale[gradstein_2020_emended$time_scale$scale==stratigraphic_scale,];
finest_chronostrat <- finest_chronostrat[finest_chronostrat$chronostratigraphic_rank==strat_rank & finest_chronostrat$interval_sr=="",];
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];
finest_chronostrat <- finest_chronostrat[finest_chronostrat$ma_lb>=min(abs(fuzzy_ranges$la_ub)),];

st <- sum(max(abs(branch_ranges$start)) <= abs(finest_chronostrat$ma_lb));	# starts!
en <- 1 + sum(min(abs(fuzzy_ranges$la_ub)) < abs(finest_chronostrat$ma_ub));	# endsd!

standard_time_scale <- finest_chronostrat[st:en,];
standard_time_scale <- standard_time_scale[!is.na(standard_time_scale$record_no),];
standard_time_scale$ma_lb <- -abs(standard_time_scale$ma_lb);
standard_time_scale$ma_ub <- -abs(standard_time_scale$ma_ub);
standard_time_scale <- standard_time_scale[standard_time_scale$interval_sr=="",];
#standard_time_scale[standard_time_scale$interval]
standard_time_scale$st[standard_time_scale$interval=="Latest Famennian"] <- "Fm4"
# Add names/symbols for plotting
strat_names <- standard_time_scale$st;
# Get the colors for time units
strat_colors <- standard_time_scale$color;
# Set the oldest and youngest intervals by names on whatever chronostratigraphic scale you used in the analysis
oldest_interval <- standard_time_scale$interval[1];
youngest_interval <- standard_time_scale$interval[nrow(standard_time_scale)];
# get the time scale that will be plotted
time_scale_to_plot <- sort(unique(c(standard_time_scale$ma_lb,standard_time_scale$ma_ub)));
time_scale_to_plot[1] <- ceiling(max(abs(branch_ranges$start)));
time_scale_to_plot[length(time_scale_to_plot)] <- -min(fuzzy_ranges$la_ub);
time_scale_to_plot <- -abs(time_scale_to_plot);
# here, I want to adjust the earliest date: I don't want the whole Cambrian
yearbreaks <- as.numeric(set_axis_breaks_new(max(abs(branch_ranges$start))-min(fuzzy_ranges$la_ub)))
#time_scale_to_plot[1] <- -ceiling(max(abs(branch_ranges$start))/max(yearbreaks))*max(yearbreaks);
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
ysize <- 10;
xsize <- 1.25*5;
mxy <- max(phylo_axis)
mny <- 0;
#mny <- min(phylo_axis);									# set maximum y value
print_names <- TRUE;
print_pix <- FALSE;
strat_label_size <- 2;					# size of the labels for chronostratigraphic units
myr_size <- 1;
#phylo_axis <- phylo_axis+(120-max(phylo_axis));
y_adj <- -max(phylo_axis)/20;
mxy <- max(phylo_axis)+y_adj;
ysize <- 7.5;
Phanerozoic_Timescale_Plot_Flexible(onset,end=end+end_extra,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma    ",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=strat_label_size,font=franky,xlab_size=1.5,myr_size=myr_size);
#### Plot Tree ####
sp <- 0;
while (sp < notu) {
#for (sp in 1:notu)	{
	sp <- sp+1;
	if (!sp %in% anagenetic_fuzzies$anc_no & !sp %in% anagenetic_fuzzies$dsc_no)	{
		# upper limit exceeds lower limit, but no definite range
		if (fuzzy_ranges$fa_lb[sp]==fuzzy_ranges$la_lb[sp] && fuzzy_ranges$fa_ub[sp]>fuzzy_ranges$la_ub[sp])	{
			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]+y_adj-1/3,-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]+y_adj+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]+y_adj-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+y_adj+1/3,col=makeTransparent(otu_colors[sp],alpha=50),border=makeTransparent(otu_colors[sp],alpha=50),lwd=0);
#			print(sp);
			} else if (fuzzy_ranges$fa_ub[sp]==fuzzy_ranges$la_ub[sp] && fuzzy_ranges$fa_lb[sp]>fuzzy_ranges$la_lb[sp])	{
		# lower limit exceeds exceeds limit, but no definite range
			rect(-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+y_adj-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+y_adj+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]+y_adj-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+y_adj+1/3,col=makeTransparent(otu_colors[sp],alpha=50),border=makeTransparent(otu_colors[sp],alpha=50),lwd=0);
#			print(sp);
			} else	{
			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]+y_adj-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+y_adj+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			}
		if (abs(fuzzy_ranges$fa_ub[sp])>abs(fuzzy_ranges$la_lb[sp]))
			rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]+y_adj-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+y_adj+1/3,col=otu_colors[sp],border=otu_colors[sp],lwd=0);
#		if (print_names)	text(-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp],fuzzy_ranges$taxon[sp],family=franky,pos=4,cex=0.5);
#		if ((fuzzy_ranges$fa_lb[sp]==fuzzy_ranges$la_lb[sp]) || (fuzzy_ranges$fa_ub[sp]==fuzzy_ranges$la_ub[sp]))	{
#			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100));
#			} else	{
#			rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100));
#			rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+1/3,col=otu_colors[sp],border=otu_colors[sp]);
#			}
		} else if (sp %in% anagenetic_fuzzies$dsc_no)	{
		fuzzy_dude <- match(sp,anagenetic_fuzzies$dsc_no);
		rect(-abs(anagenetic_fuzzies$fa_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+y_adj-1/3,-abs(anagenetic_fuzzies$la_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+y_adj+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
		if (abs(anagenetic_fuzzies$fa_ub[fuzzy_dude])>abs(anagenetic_fuzzies$la_lb[fuzzy_dude]))
			rect(-abs(anagenetic_fuzzies$fa_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+y_adj-1/3,-abs(anagenetic_fuzzies$la_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+y_adj+1/3,col=otu_colors[sp],lwd=0);
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
		segments(branch_ranges$start[nn],phylo_axis[nn]+y_adj,branch_ranges$finish[nn],phylo_axis[nn]+y_adj,lwd=1,lty=2);
		segments(branch_ranges$start[nn],phylo_axis[nn]+y_adj,branch_ranges$finish[nn],phylo_axis[nn]+y_adj,lwd=0.5,col=otu_colors[nn]);
		} else if (nn>(notu+1))	{
		nd <- nn-notu;
		if (observed_nodes_all[nd]==0)	{
			segments(branch_ranges$start[nn],phylo_axis[nn]+y_adj,branch_ranges$finish[nn],phylo_axis[nn]+y_adj,lwd=2);
			brcol <- paste("gray",round(((1-clade_posteriors[nd])*100),0),sep="");
			segments(branch_ranges$start[nn],phylo_axis[nn]+y_adj,branch_ranges$finish[nn],phylo_axis[nn]+y_adj,lwd=1,col=brcol);
#			segments(branch_ranges$start[nn],phylo_axis[nn]+y_adj,branch_ranges$finish[nn],phylo_axis[nn]+y_adj,lwd=5*clade_posteriors[nd]);
			} else	{
			segments(branch_ranges$start[nn],phylo_axis[nn]+y_adj,branch_ranges$finish[nn],phylo_axis[nn]+y_adj,lwd=1,lty=3);
#			segments(branch_ranges$start[nn],phylo_axis[nn]+y_adj,branch_ranges$finish[nn],phylo_axis[nn]+y_adj,lwd=5*clade_posteriors[nd],lty=3);
			}
		} # end case for clade
	}
for (nd in 1:nNodes)	{
	htu <- nd+notu;
	f1 <- mat_tree[nd,mat_tree[nd,]>0];
	segments(branch_ranges$finish[htu],min(phylo_axis[f1]+y_adj),branch_ranges$finish[htu],max(phylo_axis[f1]+y_adj));
#	segments(branch_ranges$finish[htu],branch_ranges$start[htu]);
	}
# plot inferred origin times of otus
for (sp in 1:notu)	points(-abs(branch_ranges$finish[sp]),phylo_axis[sp]+y_adj,pch=8,cex=10/(notu/2));
for (sp in 1:notu)	points(-abs(branch_ranges$finish[sp]),phylo_axis[sp]+y_adj,pch=8,cex=10/(notu/2),lwd=0.25,col=otu_colors[sp]);
top <- mxy+1;
otu_font_size <- max(0.25,ysize*39/(6*notu));
adj_name <- -(end+end_extra-onset)/80;
if (print_names)	{
#	otu_names <- gsub("Crassatellites","Crassatella",otu_names);
	sp <- 0;
#	for (sp in 1:notu)	{
	while (sp < notu)	{
		sp <- sp+1;
		if (!sp %in% anagenetic_fuzzies$anc_no && !sp %in% anagenetic_fuzzies$dsc_no)	{
			xx <- max(branch_ranges$finish[sp],-abs(fuzzy_ranges$la_ub[sp]));
			text(xx+adj_name,phylo_axis[sp]+y_adj,otu_names[sp],pos=4,cex=otu_font_size,family=franky);
			} else	if (sp %in% anagenetic_fuzzies$dsc_no)	{
			fuzzy_dude <- match(sp,anagenetic_fuzzies$dsc_no);
			if (is.na(fuzzy_dude))
			  fuzzy_dude <- match(sp,anagenetic_fuzzies$anc_no);
#		lineage_name <- paste(anagenetic_fuzzies$ancestor[fuzzy_dude],"->",anagenetic_fuzzies$descendant[fuzzy_dude],sep="")
			xx <- max(branch_ranges$finish[sp],-abs(fuzzy_ranges$la_ub[sp]));
			text(xx+adj_name,phylo_axis[sp]+y_adj,anagenetic_fuzzies$taxon[fuzzy_dude],pos=4,cex=otu_font_size,family=franky);
			}
		}
	}

# print key: ####
#ysize <- 3;
y <- 3;
Phanerozoic_Timescale_Plot_Flexible(onset,end=end+end_extra,time_scale_to_plot,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names,strat_colors,plot_title,ordinate,abscissa="Ma    ",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size=2*strat_label_size,font=franky,xlab_size=1.5,myr_size=myr_size);
ngroups <- length(group_colors)
yspan <- mxy-mny;
xspan <- end-onset;
py <- mxy;
dy <- 1.5*yspan/40;

for (i in 1:ngroups)	{
	if (variable_to_plot=="taxonomy") {
		j <- match(ingroup_ht[i],names(group_colors));
		} else if (variable_to_plot=="geography")	{
		j <- i;
		}
	rect(onset,py,onset+xspan/30,py-yspan/40,col=makeTransparent(group_colors[j],alpha=100),border=makeTransparent(group_colors[j],alpha=100));
	rect(onset+xspan/30,py,onset+xspan/15,py-yspan/40,col=group_colors[j],border=group_colors[j]);
	text(onset+xspan/15,py-yspan/80,paste(":",names(group_colors)[j]),family=franky,pos=4);
	py <- py-dy;
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
#	rect(hpd$lb[cfn],phylo_axis[cfn]-0.125,hpd$ub[cfn],phylo_axis[cfn]+0.125,col="white",lwd=0);
#	if (i==1)	{
#		temp_col <- makeTransparent("magenta",alpha=50);
#		rect(hpd$lb[cfn],phylo_axis[cfn]-0.125,hpd$ub[cfn],phylo_axis[cfn]+0.125,col=temp_col,lwd=0);
#		} else	{
		temp_col <- makeTransparent(higher_taxon_colors[i],alpha=100);
		rect(hpd$lb[cfn],phylo_axis[cfn]-0.125,hpd$ub[cfn],phylo_axis[cfn]+0.125,col=temp_col,lwd=0);
#		}
	}
#for (pc in 1:length(pictured))	{
#	jpg_file <- paste(this_directory,pictured[pc],"_horizontal.jpg",sep="");
#	txn <- match(pictured[pc],otu_names_nex);
#	add_jpeg(p,x=-abs(fuzzy_ranges$la_ub[txn])+x_add,y=phylo_axis[txn]-0.5,add=TRUE)
#	}

{
mctree2 <- scan(file=mctree_file2,what=character(),sep="\n")[1];
mctree2_info <- read_newick_string_with_branch_lengths(mctree2);
branch_durations2 <- mctree2_info$branch_lengths; # sum(is.na(branch_durations2))
vector_tree2 <- mctree2_info$vector_tree; # sum(is.na(vector_tree2))
htu_heights2 <- mctree2_info$htu_heights; # sum(is.na(htu_heights2))
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

taxlabels <- match("taxlabels",tolower(mctree));
otu_names_nex <- mctree[(taxlabels+1):length(mctree)];
otu_names_nex <- otu_names_nex[1:(match(";",otu_names_nex)-1)];
otu_names <- gsub("_"," ",otu_names_nex);
#print(otu_names[2]);

#print(otu_names[2]);
	
notu <- length(otu_names);
#print(otu_names[2]);


# make sure that ranges & taxa are in the same order!
#fuzzy_ranges <- fuzzy_ranges[match(otu_names,fuzzy_ranges$taxon),];
#otu_information <- otu_information[match(otu_names,otu_information$taxon),];

tree_line <- 1+match("begin trees;",tolower(mctree));
tree_line <- 1+match("translate;",tolower(mctree));
tree_line <- (1:length(mctree))[gsub("TREE","",mctree)!=mctree];
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&R\\] ","",mctree[tree_line]);
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&R\\]","",mctree[tree_line]);
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&U\\] ","",newick_string_full);
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&U\\]","",newick_string_full);
newick_string_full <- tree_info <- gsub("tree TREE1 = ","",newick_string_full);
newick_string_full <- tree_info <- gsub("tree TREE = ","",newick_string_full);
newick_string_full <- tree_info <- gsub(" \\(","\\(",newick_string_full);
if (sum(sapply(otu_names,is.subspecies))>0)	{
	subspecies <- otu_names[(1:notu)[sapply(otu_names,is.subspecies)]];
	new_species <- sapply(subspecies,elevate_subspecies_to_species);
	new_species_nex <- gsub(" ","_",new_species);
	subspecies_nex <- gsub(" ","_",subspecies);
	otu_names_nex[otu_names_nex %in% subspecies_nex] <- new_species_nex;
	otu_names_orig <- otu_names;
	otu_names[otu_names %in% subspecies] <- new_species;
	} else	{
	otu_names_orig <- otu_names;
	}
otu_names_orig_nex <- gsub(" ","_",otu_names_orig);
#otu_names[2]
if (gsub(otu_names_nex[1],"",newick_string_full)==newick_string_full)
	for (i in notu:1)
		newick_string_full <- gsub(paste(i,"\\[",sep=""),paste(otu_names_nex[i],"\\[",sep=""),newick_string_full);

mctree_info <- read_newick_string_mcmc_beast(newick_string_full,otu_names);
clade_posteriors1 <- mctree_info$clade_posteriors;
names(mctree_info$ancestral)[!names(mctree_info$ancestral) %in% otu_names_orig_nex] <- otu_names_orig_nex[!names(mctree_info$ancestral) %in% otu_names_orig_nex];
mctree_info$ancestral <- mctree_info$ancestral[!is.na(names(mctree_info$ancestral))];
names(mctree_info$branch_durations)[1:notu] <- names(mctree_info$ancestral);

vector_tree1a <- mctree_info$vector_tree;
branch_durations <- branch_durations1 <- mctree_info$branch_durations;
names(vector_tree1a) <- names(branch_durations);
clade_posteriors2 <- clade_posteriors <- clade_posteriors1;
clade_posteriors2[clade_posteriors2>0] <- 0;
vector_tree1 <- vector_tree1a[order(names(vector_tree1a[1:notu]))];
names(vector_tree1) <- names(vector_tree1a[1:notu][order(names(vector_tree1a[1:notu]))]);
vector_tree1a <- vector_tree1a[(notu+1):length(vector_tree1a)];
vector_tree <- vector_tree1 <- c(vector_tree1,vector_tree1a);
for (i in 1:length(vector_tree))	{
	if (vector_tree1[i]!=vector_tree2[i])	{
		nd1 <- vector_tree1[i];
		nd2 <- vector_tree2[i];
		clade_posteriors2[nd2-notu] <- clade_posteriors1[nd1-notu];
		vector_tree[i] <- nd2;
		}
	}


}