#### Setup Program ####
this_directory <- "~/Documents/R_Projects/Tree_Drawing/";
setwd(this_directory);
common_source_folder <- "~/Documents/R_Projects/Common_R_Source_Files/";	# directory to folder where you keep common source
franky <- "Franklin Gothic Medium";  #Or whatever other R-compatible font you'd like to use...

#install.packages("png")
library(png);
library(zoo);
source(paste(common_source_folder,"General_Plot_Templates.r",sep=""));
source(paste(common_source_folder,"Nexus_File_Routines.r",sep=""));

data_for_R_folder <- "~/Documents/R_Projects/Data_for_R/";
load(paste(data_for_R_folder,"Gradstein_2020_Augmented.RData",sep="")); # refined Gradstein 2012 timescale & biozonations
gradstein_2020_emended$time_scale$scale[gradstein_2020_emended$time_scale$interval=="Early Cambrian"] <- "International Old";

merge_anagenetics <- FALSE;
tree_to_print <- "3U";
# set parameters & files for printing tree ####
if (tree_to_print %in% c("1L","1U"))	{
	strat_rank <- "Epoch";		# I recommend either "Stage" or "Stage Slice"
	rank_to_plot <- "order";	# choose which taxonomic rank you want to print; irrelevant if you don't choose plot_taxonomy <- T;
	range_file <- paste(this_directory,"taxon_info/Diploporita_Fuzzy_Ranges.xlsx",sep="");
	otu_information_file <- paste(this_directory,"taxon_info/Diploporita_Information.xlsx",sep="");
	if (tree_to_print=="1L")	{
		mctree_file <- paste(tree_drawing_folder,"trees/Diploplorites_no_outgroup.tre",sep="");
#		mctree_file <- paste(this_directory,"trees/dip_mcc.tre",sep="");
		} else if (tree_to_print=="1U")	{
#		mctree_file <- paste(this_directory,"trees/diploploranNoOutGroupNoJump.mcc.tre",sep="");
		mctree_file <- paste(tree_drawing_folder,"trees/dip_noOG.mcc.tre",sep="");
		mctree_file <- paste(tree_drawing_folder,"trees/diploporanNoOutGroupNoJump.mcc.tre",sep="");
		}
	} else if (tree_to_print %in% c("2U","2L"))	{
	strat_rank <- "Epoch";		# I recommend either "Stage" or "Stage Slice"
	rank_to_plot <- "order";	# choose which taxonomic rank you want to print; irrelevant if you don't choose plot_taxonomy <- T;
	range_file <- paste(this_directory,"taxon_info/Blastoidea_Fuzzy_Ranges.xlsx",sep="");
	otu_information_file <- paste(this_directory,"taxon_info/Blastoidea_Information.xlsx",sep="");
	if (tree_to_print=="2U")	{
		mctree_file <- paste(this_directory,"trees/blastoid_constrained.mcc.tre",sep="");
		mctree_file <- paste(tree_drawing_folder,"trees/blastoids_clade_constraints_unlinked.tre",sep="");
		} else if (tree_to_print=="2L")	{
	#	mctree_file <- paste(tree_drawing_folder,"trees/blastoid.tre",sep="");
		mctree_file <- paste(tree_drawing_folder,"trees/blastoids_clade_constraints_linked.tre",sep="");
#		mctree_file <- paste(tree_drawing_folder,"Blastoids_clade_constraints.tre",sep="");
		}
	} else if (tree_to_print %in% c("3U","3L"))	{
	range_file <- paste(tree_drawing_folder,"Paracrinoidea_Fuzzy_Ranges.xlsx",sep="");
	taxon_information_file <- paste(tree_drawing_folder,"Paracrinoidea_Taxon_Information.xlsx",sep="");
#	mctree_file <- paste(tree_drawing_folder,"paracrinoid_unconstrained_mcc.tre",sep="");
	strat_rank <- "Stage";		# I recommend either "Stage" or "Substage"
	rank_to_plot <- "family";	# choose which taxonomic rank you want to print; irrelevant if you don't choose plot_taxonomy <- TRUE;
	if (tree_to_print=="3U")	{
		mctree_file <- paste(tree_drawing_folder,"trees/ParacrinoidsEstonocystisUnlinked.tre",sep="");
		} else if (tree_to_print=="3L")	{
		mctree_file <- paste(tree_drawing_folder,"trees/NoCheiro_mcc-2024-07-30.tre",sep="");
		mctree_file <- paste(tree_drawing_folder,"trees/CheiroOutgroup_mcc-2024-07-30.tre",sep="");
		}
	rank_to_plot <- "family";	# choose which taxonomic rank you want to print; irrelevant if you don't choose plot_taxonomy <- T;
	range_file <- paste(this_directory,"taxon_info/Paracrinoidea_Fuzzy_Ranges.xlsx",sep="");
	otu_information_file <- paste(this_directory,"taxon_info/Paracrinoidea_Information.xlsx",sep="");
	}

# load stratigraphic data & set it up for figures ####
if (gsub(".csv","",otu_information_file)!=otu_information_file)	{
	otu_information <- read.csv(otu_information_file,header=T,stringsAsFactors = F);
	} else if (gsub(".xlsx","",otu_information_file)!=otu_information_file)	{
	otu_information <- as.data.frame(readxl::read_xlsx(otu_information_file));
	}
colnames(otu_information)[colnames(otu_information)=="species"] <- "taxon";
colnames(otu_information) <- tolower(colnames(otu_information));

if (gsub(".csv","",range_file)!=range_file)	{
	fuzzy_ranges <- read.csv(range_file,header=T,stringsAsFactors = F);
	} else if (gsub(".xlsx","",range_file)!=range_file)	{
	fuzzy_ranges <- as.data.frame(readxl::read_xlsx(range_file));
	}
colnames(fuzzy_ranges)[colnames(fuzzy_ranges)=="species"] <- "taxon";
colnames(fuzzy_ranges) <- tolower(colnames(fuzzy_ranges));
fuzzy_ranges$taxon <-gsub("_"," ",fuzzy_ranges$taxon);
#fuzzy_ranges$taxon <- gsub("Crassatellites","Crassatella",fuzzy_ranges$taxon);
faux_recent <- min(fuzzy_ranges$fa_ub);
oldest_fa <- -(max(fuzzy_ranges$fa_ub)+max(fuzzy_ranges$fa_lb))/2;
rescaled_fa_bounds <- fuzzy_ranges;
rescaled_fa_bounds$la_lb <- rescaled_fa_bounds$la_ub <- NULL;
rescaled_fa_bounds$fa_lb <- rescaled_fa_bounds$fa_lb - faux_recent;
rescaled_fa_bounds$fa_ub <- rescaled_fa_bounds$fa_ub - faux_recent;

finest_chronostrat <- subset(gradstein_2020_emended$time_scale,gradstein_2020_emended$time_scale$chronostratigraphic_rank==strat_rank)
finest_chronostrat <- subset(finest_chronostrat,finest_chronostrat$scale=="International")
finest_chronostrat <- finest_chronostrat[finest_chronostrat$interval_sr=="",];
finest_chronostrat <- finest_chronostrat[order(-finest_chronostrat$ma_lb),];

stage_slice_chronostrat <- subset(gradstein_2020_emended$time_scale,gradstein_2020_emended$time_scale$scale=="Stage Slice");
stage_slice_chronostrat <- stage_slice_chronostrat[order(-abs(stage_slice_chronostrat$ma_lb)),];

if (strat_rank=="Stage Slice")	finest_chronostrat <- stage_slice_chronostrat;

#### Get Tree Information for Plotting ####
mctree <- scan(file=mctree_file,what=character(),sep="\n");
mctree <- gsub("\t","",mctree)
begin_taxa <- match("begin taxa;",tolower(mctree));
notu <- as.numeric(strsplit(gsub(";","",mctree[begin_taxa+1]),split = "=")[[1]][2]);

taxlabels <- match("taxlabels",tolower(mctree));

fuzzy_ranges <- fuzzy_ranges[fuzzy_ranges$taxon %in% gsub("_"," ",mctree[(taxlabels+1):(taxlabels+notu)]),]
#otu_names <- fuzzy_ranges$taxon[fuzzy_ranges$taxon %in% gsub("_"," ",mctree[(taxlabels+1):(taxlabels+notu)])];
otu_names <- fuzzy_ranges$taxon;
otu_names_nex <- gsub(" ","_",otu_names);
tree_line <- 1+match("begin trees;",tolower(mctree));
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&R\\]","",mctree[tree_line]);
newick_string_full <- tree_info <- gsub("tree TREE1 = \\[&U\\]","",newick_string_full);

tree_output <- read_newick_string_mcmc(newick_string_full,otu_names);
clade_posteriors <- tree_output$clade_posteriors;
sampled_ancestors <- tree_output$ancestral;

basal_taxon <- fuzzy_ranges$taxon[match(max(fuzzy_ranges$fa_lb),fuzzy_ranges$fa_lb)];
basal_date_raw <- max(((tree_output$hpd$lb+tree_output$hpd$ub)/2)[1:notu]);	# base from Bayesian tree before rescaling 0 to actual upper bound
basal_taxon_no <- match(basal_date_raw,((tree_output$hpd$lb+tree_output$hpd$ub)/2)[1:notu]);

date_rescale <- min(fuzzy_ranges$fa_ub);
fuzzy_ranges$hpd_lb <- tree_output$hpd_age$lb[1:notu];
fuzzy_ranges$hpd_ub <- tree_output$hpd_age$ub[1:notu];
#write.csv(fuzzy_ranges,"WTF.csv",row.names = FALSE);
hpd <- -(tree_output$hpd+date_rescale);

#hpd <- data.frame(-min(fuzzy_ranges$fa_lb)-tree_output$hpd);
hpd$md <- (hpd$lb+hpd$ub)/2;
branch_durations <- abs(tree_output$branch_durations);
#branch_durations[(1:notu)[sampled_ancestors==1]] <- 0;
vector_tree <- tree_output$vector_tree;
notu <- match(-1,tree_output$vector_tree)-1;
mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree=tree_output$vector_tree);
nNodes <- nrow(mat_tree);
lowest_taxon <- min(which(mat_tree %in% (1:notu),arr.ind = T));
outgp <- mat_tree[lowest_taxon,mat_tree[lowest_taxon,]<notu];
#-max(abs(hpd$md[outgp]))

### Construct proper phylogeny ####
root_ma <- min(hpd$md);
branch_ranges <- data.frame(start=as.numeric(vector_tree),
														finish=as.numeric(vector_tree));
rownames(branch_ranges) <- names(branch_durations)[!is.na(names(branch_durations))];
branch_ranges$finish <- hpd$md;

branch_ranges[notu+1,] <- root_ma;
branch_ranges$start[1:notu] <- hpd$md[vector_tree[1:notu]];
branch_ranges$start[(notu+2):(notu+nNodes)] <-hpd$md[vector_tree[(notu+2):(notu+nNodes)]];

ancestral_spc <- (1:notu)[sampled_ancestors==1];
names(ancestral_spc) <- names(sampled_ancestors)[sampled_ancestors==1];
anagenetic_ancestors <- rep(0,notu);
names(anagenetic_ancestors) <- names(sampled_ancestors);
if (length(ancestral_spc)>0)	{
	pbdb_finds <- pbdb_data_list$pbdb_finds;
	pbdb_taxonomy <- pbdb_data_list$pbdb_taxonomy;
	venn_tree <- transform_vector_tree_to_venn_tree(vector_tree);
	mat_tree <- transform_vector_tree_to_matrix_tree(vector_tree);
#	cinctan_finds <- read.csv(paste(this_directory,"Cincta_Finds.csv",sep=""),header=T,stringsAsFactors = F);
	for (an in 1:length(ancestral_spc))	{
		anc <- ancestral_spc[an];
		daughter_node <- max(which(venn_tree==anc,arr.ind = T)[,1]);
		its_node <- which(mat_tree==anc,arr.ind = T)[,1]
		desc <- mat_tree[its_node,mat_tree[its_node,]!=anc];
	#	hpd[anc,]
	#	hpd[desc,]
	#	fuzzy_ranges[anc,]
	#	fuzzy_ranges[desc,]
		# prior routine using venn_tree
		progeny <- venn_tree[daughter_node,venn_tree[daughter_node,]>0];
		progeny <- progeny[progeny!=anc];
		anc_find_data <- subset(pbdb_finds,pbdb_finds$accepted_name==otu_names[anc])
		if(sum(abs(fuzzy_ranges$la_lb[anc])<=abs(fuzzy_ranges$fa_ub[progeny]))==0)	{
#			find_data <- subset(pbdb_finds,pbdb_finds$accepted_name %in% otu_names[c(anc,progeny)]);
			accepted_names <- unique(pbdb_taxonomy$accepted_name[pbdb_taxonomy$taxon_name %in% otu_names[c(anc,progeny)]]);
			accepted_names <- accepted_names[!is.na(accepted_names)];
			poss_names <- pbdb_taxonomy$taxon_name[pbdb_taxonomy$accepted_name %in% accepted_names];
			poss_taxon_nos <- pbdb_taxonomy$taxon_no[pbdb_taxonomy$accepted_name %in% accepted_names];
			find_data <- pbdb_finds[pbdb_finds$accepted_no %in% poss_taxon_nos,]
			find_data$accepted_name <- pbdb_taxonomy$accepted_name[match(find_data$accepted_no,pbdb_taxonomy$taxon_no)];
#			pbdb_finds[pbdb_finds$identified_name=="Crassatella tumidula",]
			coccur_mat <- accersi_cooccurence_matrix(find_data);
			if (sum(coccur_mat[match(otu_names[anc],rownames(coccur_mat)),])==1)
#			if (coccur_mat[1,2]==0)
				anagenetic_ancestors[anc] <- 1;
			}
		}
	}

#nNodes <- nrow(mat_tree)
phylo_axis <- get_phylo_axis_from_newick_string_w_anagenesis(newick_string = tree_output$newick,sampled_ancestors=sampled_ancestors,anagenetic_ancestors=anagenetic_ancestors);
write.csv(cbind(sampled_ancestors,anagenetic_ancestors),"Ancestors2.csv",row.names = TRUE)
names(phylo_axis) <- c((otu_names),paste("nd_",1:nNodes,sep=""));


### use this to identify which sequences get lumped ####
# get anagenetic series of species
anagenetic_anc_no <- (1:notu)[anagenetic_ancestors==1];
anagenetic_ancestral_taxa <- names(phylo_axis[1:notu])[anagenetic_anc_no];
names(anagenetic_anc_no) <- names(anagenetic_ancestors)[anagenetic_anc_no];
axis_ranks <- rank(phylo_axis[1:notu]);
axis_ranks <- axis_ranks-min(axis_ranks)+1;
overlaps <- hist(match(axis_ranks,unique(axis_ranks)),breaks=0:length(unique(axis_ranks)),plot=FALSE)$counts
anagenetic_overlaps <- hist(phylo_axis[1:notu],breaks=(min(phylo_axis[1:notu])-1):ceiling(max(phylo_axis[1:notu])),plot=F)$counts;
phylo_axis_ana <- (1:length(anagenetic_overlaps))[anagenetic_overlaps>1] - (1-min(phylo_axis[1:notu])); # NOTE: THIS HAS TO GO DOWN TO ZERO!!!!

anagenetic_fuzzies <- data.frame(ancestor=as.character(),descendant=as.character(),
								 fa_lb=as.numeric(),fa_ub=as.numeric(),la_lb=as.numeric(),la_ub=as.numeric(),
								 anc_no=as.numeric(),dsc_no=as.numeric(),spc_sequence=as.character(),
								 phylo_axis=as.numeric(),taxon=as.character(),stringsAsFactors = F);
af <- 0;
while (af < length(phylo_axis_ana))	{
	af <- af+1;
	dummy <- data.frame(ancestor=as.character(""),descendant=as.character(""),
											fa_lb=as.numeric(0),fa_ub=as.numeric(0),la_lb=as.numeric(0),la_ub=as.numeric(0),
											anc_no=as.numeric(0),dsc_no=as.numeric(0),spc_sequence=as.character(""),
											phylo_axis=as.numeric(0),taxon=as.character(""),stringsAsFactors = F);
	ana_spc_nos <- (1:notu)[phylo_axis[1:notu] %in% phylo_axis_ana[af]];
	this_lineage <- names(phylo_axis[1:notu])[phylo_axis[1:notu] %in% phylo_axis_ana[af]];
#	this_lineage <- fuzzy_ranges$taxon[ana_spc_nos];
	dummy$ancestor <- this_lineage[this_lineage %in% anagenetic_ancestral_taxa];
	dummy$ancestor <- gsub("_"," ",dummy$ancestor);
	dummy$fa_lb <- fuzzy_ranges$fa_lb[fuzzy_ranges$taxon==dummy$ancestor]
	dummy$fa_ub <- fuzzy_ranges$fa_ub[fuzzy_ranges$taxon==dummy$ancestor]
	descendants <- gsub("_"," ",this_lineage[!this_lineage %in% anagenetic_ancestral_taxa]);
	dummy$la_lb <- min(fuzzy_ranges$la_lb[match(descendants,fuzzy_ranges$taxon)])
	dummy$la_ub <- min(fuzzy_ranges$la_ub[match(descendants,fuzzy_ranges$taxon)])
	dummy$descendant <- gsub("_"," ",paste(this_lineage[!this_lineage %in% anagenetic_ancestral_taxa],collapse=", "));
	dummy$phylo_axis <- phylo_axis_ana[af];
	dummy$anc_no <- match(dummy$ancestor,fuzzy_ranges$taxon);
	dummy$dsc_no <- match(dummy$descendant,fuzzy_ranges$taxon);
	this_lineage <- gsub("_"," ",this_lineage);
	this_lineage <- c(dummy$ancestor,this_lineage[!this_lineage %in% dummy$ancestor]);
	dummy$spc_sequence <- paste(ana_spc_nos[order(branch_ranges$finish[ana_spc_nos])],collapse=",");
	dummy$taxon <- paste(this_lineage,collapse="->");
	anagenetic_fuzzies <- rbind(anagenetic_fuzzies,dummy);
	}

#### plot time scale ####
basal_divergence <- ceiling(max(abs(branch_ranges$start))/5)*5;
if (tree_to_print=="1U")	basal_divergence <- 540;
if (tree_to_print=="2U")	basal_divergence <- 503.337;
if (tree_to_print=="3U")	basal_divergence <- 560.0;
st <- sum(basal_divergence <= abs(finest_chronostrat$ma_lb));	# starts!
en <- sum(min(abs(fuzzy_ranges$la_ub)) <= abs(finest_chronostrat$ma_lb));	# endsd!
if (tree_to_print=="1L")	en <- 14;

standard_time_scale <- finest_chronostrat[st:en,];
standard_time_scale$ma_lb <- -abs(standard_time_scale$ma_lb);
standard_time_scale$ma_ub <- -abs(standard_time_scale$ma_ub);
standard_time_scale <- standard_time_scale[standard_time_scale$interval_sr=="",]
standard_time_scale$st[standard_time_scale$interval %in% "Avalon"] <- "Ed3";
standard_time_scale$st[standard_time_scale$interval %in% "White Sea"] <- "Ed4";
standard_time_scale$st[standard_time_scale$interval %in% "Nama"] <- "Ed5";

time_scale <- gradstein_2020_emended$time_scale[gradstein_2020_emended$time_scale$chronostratigraphic_rank %in% "Period" & gradstein_2020_emended$time_scale$scale %in% "International",];
time_scale$ma_lb <- -time_scale$ma_lb; time_scale$ma_ub <- -time_scale$ma_ub;
time_scale <- time_scale[order(time_scale$ma_lb),];
time_scale <- time_scale[time_scale$ma_ub>min(standard_time_scale$ma_lb),];
time_scale <- time_scale[time_scale$ma_lb<max(standard_time_scale$ma_ub),];

# Add names/symbols for plotting
strat_names_minor <- standard_time_scale$st;
strat_names_major <- time_scale$interval;
# Get the colors for time units
strat_colors_minor <- standard_time_scale$color;
names(strat_colors_minor) <- strat_names_minor;
strat_colors_major <- time_scale$color;
# Set the oldest and youngest intervals by names on whatever chronostratigraphic scale you used in the analysis
oldest_interval <- standard_time_scale$interval[1];
youngest_interval <- standard_time_scale$interval[nrow(standard_time_scale)];
# get the time scale that will be plotted
time_scale_to_plot_minor <- unique(c(standard_time_scale$ma_lb,standard_time_scale$ma_ub));
# here, I want to adjust the earliest date: I don't want the whole Cambrian
#time_scale_to_plot_minor[1] <- -abs(stage_slice_chronostrat$ma_lb[sum(-abs(stage_slice_chronostrat$ma_lb) < (-ceiling(max(abs(branch_ranges$start)))))]);
time_scale_to_plot_minor[1] <- -abs(stage_slice_chronostrat$ma_lb[sum(-abs(stage_slice_chronostrat$ma_lb) < (-ceiling(abs(basal_divergence))))]);
time_scale_to_plot_minor[1] <- -ceiling(abs(basal_divergence));
time_scale_to_plot_minor[length(time_scale_to_plot_minor)] <- -floor(min(abs(fuzzy_ranges$la_ub))/5)*5;

time_scale_to_plot_major <- unique(c(time_scale$ma_lb,time_scale$ma_ub));
time_scale_to_plot_major[1] <- time_scale_to_plot_minor[1];
time_scale_to_plot_major[time_scale_to_plot_major==max(time_scale_to_plot_major)] <- max(time_scale_to_plot_minor);

# Now, set the onset and end of the x-axis
onset <- min(time_scale_to_plot_minor);
end <- max(time_scale_to_plot_minor);
# Finally, set up breaks for x-axis
#yearbreaks <- c(5,25,50);					# set breaks for x-axis (minor, medium & major)
yearbreaks <- as.numeric(set_axis_breaks(max_no=end,min_no=onset));
# now, set up the y-axis: this will reflect your data

use_strat_labels <- T;						# if T, then strat_names_minor will be plotted on X-axis inside boxes
alt_back <- F;								# if T, then the background will alternat shades between major intervals
plot_title <- "";							# Name of the plot; enter "" for nothing
ordinate <- "";								# Label of Y-axis
hues <- "T";								# If T, then IGN stratigraphic colors will be used
colored <- "base";							# Where IGN stratigraphic colors should go
mxy <- max(phylo_axis)
mny <- min(phylo_axis);									# set maximum y value
#mny <- 0;
yearbreaks <- sort(as.numeric(set_axis_breaks_new(max(abs(time_scale_to_plot_minor))-min(abs(time_scale_to_plot_minor)))));
#yearbreaks <- c(1,0.5,0.1);
yearbreaks <- sort(yearbreaks);

end_extra <- abs(end-onset)/4;
ysize <-  min(6,4.285714285*notu/20);
#ysize <- 5;
xsize <- 1.25*5; 
mxy <- max(phylo_axis)
mny <- min(phylo_axis);									# set maximum y value
ysize <- xsize/1.25;
print_names <- TRUE;
print_pix <- FALSE;
strat_label_size <- 2/3;					# size of the labels for chronostratigraphic units
myr_size <- 1;
strat_label_size_minor <- 1.5;
strat_label_size_major <- 0.85;
if (onset > -510)	strat_names_major[strat_names_major=="Cambrian"] <- "Cam";
if (onset > -570)	strat_names_major[strat_names_major=="Ediacaran"] <- "Ediac";
strat_names_minor[strat_names_minor=="lEd"] <- "L.Ed";
strat_names_minor[strat_names_minor=="Dpn"] <- "Dp";
#strat_names_minor[strat_names_minor %in% c("Dp","Dpn")] <- "";
if (end < -270)	strat_names_major[strat_names_major=="Permian"] <- "Perm";
Phanerozoic_Timescale_Plot_Hierarchical(onset,end+end_extra,time_scale_to_plot_minor,time_scale_to_plot_major,mxy,mny=(mny-(mxy-mny)/20),use_strat_labels,strat_names_minor,strat_names_major,strat_colors_minor,strat_colors_major,plot_title,ordinate,abscissa="Ma                   ",yearbreaks,xsize=xsize,ysize=ysize,hues=hues,colored=colored,alt_back=alt_back,strat_label_size_minor=strat_label_size_minor,strat_label_size_major=strat_label_size_major,xlab_size=1.5,myr_size=myr_size,font=franky);
#### Plot Tree ####
otu_colors <- otu_information$color;
names(otu_colors) <- otu_information$taxon;
# plot stratigraphic ranges
if (merge_anagenetics)	{
	for (sp in 1:notu)	{
	#	sp <- sp+1;
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
#			if ((fuzzy_ranges$fa_lb[sp]==fuzzy_ranges$la_lb[sp]) || (fuzzy_ranges$fa_ub[sp]==fuzzy_ranges$la_ub[sp]))	{
#				rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100));
#				} else	{
#				rect(-abs(fuzzy_ranges$fa_lb[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_ub[sp]),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100));
#				rect(-abs(fuzzy_ranges$fa_ub[sp]),phylo_axis[sp]-1/3,-abs(fuzzy_ranges$la_lb[sp]),phylo_axis[sp]+1/3,col=otu_colors[sp],border=otu_colors[sp]);
#				}
			} else if (sp %in% anagenetic_fuzzies$dsc_no)	{
			fuzzy_dude <- match(sp,anagenetic_fuzzies$dsc_no);
			rect(-abs(anagenetic_fuzzies$fa_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]-1/3,-abs(anagenetic_fuzzies$la_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
			if (abs(anagenetic_fuzzies$fa_ub[fuzzy_dude])>abs(anagenetic_fuzzies$la_lb[fuzzy_dude]))
				rect(-abs(anagenetic_fuzzies$fa_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]-1/3,-abs(anagenetic_fuzzies$la_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+1/3,col=otu_colors[sp],lwd=0);
#		if ()
			}
		}
	} else	{
	for (sp in 1:notu)	{
	#	sp <- sp+1;
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
			} else if (sp %in% anagenetic_fuzzies$anc_no)	{
			fuzzy_dude <- match(sp,anagenetic_fuzzies$anc_no);
			sep_species <- as.numeric(strsplit(anagenetic_fuzzies$spc_sequence[fuzzy_dude],",")[[1]]);
			relv_ranges <- fuzzy_ranges[sep_species,];
			overlapping_ranges <- c();
			overlap_matrix <- array(0,dim=c(length(sep_species),length(sep_species)));
			colnames(overlap_matrix) <- rownames(overlap_matrix) <- sep_species;
			for (i in 1:(nrow(relv_ranges)-1))	{
				overlap_matrix[i,i] <- 1;
				for (j in (i+1):nrow(relv_ranges))	{
					overlap_matrix[j,j] <- 1;
					if (relv_ranges$la_ub[i]<relv_ranges$fa_lb[j])	{
						overlap_matrix[i,j] <- overlap_matrix[j,i] <- 1;
						overlapping_ranges <- unique(c(overlapping_ranges,c(sep_species[i],sep_species[j])));
						}
					}
				}
			disjunct_ranges <- sep_species[!sep_species %in% overlapping_ranges];
			for (ss in 1:length(sep_species))	{
				if (sep_species[ss] %in% disjunct_ranges)	{
					spp <- sep_species[ss];
					if (fuzzy_ranges$fa_lb[spp]==fuzzy_ranges$la_lb[spp] && fuzzy_ranges$fa_ub[spp]>fuzzy_ranges$la_ub[spp])	{
						rect(-abs(fuzzy_ranges$fa_lb[spp]),phylo_axis[spp]-1/3,-abs(fuzzy_ranges$fa_ub[spp]),phylo_axis[spp]+1/3,col=makeTransparent(otu_colors[spp],alpha=100),border=makeTransparent(otu_colors[spp],alpha=100),lwd=0);
						rect(-abs(fuzzy_ranges$fa_ub[spp]),phylo_axis[spp]-1/3,-abs(fuzzy_ranges$la_ub[spp]),phylo_axis[spp]+1/3,col=makeTransparent(otu_colors[spp],alpha=50),border=makeTransparent(otu_colors[spp],alpha=50),lwd=0);
						} else if (fuzzy_ranges$fa_ub[spp]==fuzzy_ranges$la_ub[spp] && fuzzy_ranges$fa_lb[spp]>fuzzy_ranges$la_lb[spp])	{
					# lower limit exceeds exceeds limit, but no definite range
						rect(-abs(fuzzy_ranges$la_lb[spp]),phylo_axis[spp]-1/3,-abs(fuzzy_ranges$la_ub[spp]),phylo_axis[spp]+1/3,col=makeTransparent(otu_colors[spp],alpha=100),border=makeTransparent(otu_colors[spp],alpha=100),lwd=0);
						rect(-abs(fuzzy_ranges$fa_lb[spp]),phylo_axis[spp]-1/3,-abs(fuzzy_ranges$la_lb[spp]),phylo_axis[spp]+1/3,col=makeTransparent(otu_colors[spp],alpha=50),border=makeTransparent(otu_colors[spp],alpha=50),lwd=0);
						} else	{
						rect(-abs(fuzzy_ranges$fa_lb[spp]),phylo_axis[spp]-1/3,-abs(fuzzy_ranges$la_ub[spp]),phylo_axis[spp]+1/3,col=makeTransparent(otu_colors[spp],alpha=100),border=makeTransparent(otu_colors[spp],alpha=100),lwd=0);
						}
					if (abs(fuzzy_ranges$fa_ub[spp])>abs(fuzzy_ranges$la_lb[spp]))
						rect(-abs(fuzzy_ranges$fa_ub[spp]),phylo_axis[spp]-1/3,-abs(fuzzy_ranges$la_lb[spp]),phylo_axis[spp]+1/3,col=otu_colors[spp],border=otu_colors[spp],lwd=0);
					} else	{
					if (sum(overlap_matrix)==(nrow(overlap_matrix)*ncol(overlap_matrix)))	{
						this_sequence <- sep_species;
						} else	{
						
						}
					this_seq_ranges <- fuzzy_ranges[this_sequence,];
					lineage_range <- this_seq_ranges[1,];
					lineage_range$la_lb <- min(lineage_range$la_lb);
					lineage_range$la_ub <- min(lineage_range$la_ub);

					if (lineage_range$fa_lb==lineage_range$la_lb && lineage_range$fa_ub>lineage_range$la_ub)	{
						rect(-abs(lineage_range$fa_lb),phylo_axis[sp]-1/3,-abs(lineage_range$fa_ub),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
						rect(-abs(lineage_range$fa_ub),phylo_axis[sp]-1/3,-abs(lineage_range$la_ub),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=50),border=makeTransparent(otu_colors[sp],alpha=50),lwd=0);
						} else if (lineage_range$fa_ub==lineage_range$la_ub && lineage_range$fa_lb>lineage_range$la_lb)	{
				# lower limit exceeds exceeds limit, but no definite range
						rect(-abs(lineage_range$la_lb),phylo_axis[sp]-1/3,-abs(lineage_range$la_ub),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
						rect(-abs(lineage_range$fa_lb),phylo_axis[sp]-1/3,-abs(lineage_range$la_lb),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=50),border=makeTransparent(otu_colors[sp],alpha=50),lwd=0);
						} else	{
						rect(-abs(lineage_range$fa_lb),phylo_axis[sp]-1/3,-abs(lineage_range$la_ub),phylo_axis[sp]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
						}
					if (abs(lineage_range$fa_ub)>abs(lineage_range$la_lb))
						rect(-abs(lineage_range$fa_ub),phylo_axis[sp]-1/3,-abs(lineage_range$la_lb),phylo_axis[sp]+1/3,col=otu_colors[sp],border=otu_colors[sp],lwd=0);
					}
				
				}
#			rect(-abs(anagenetic_fuzzies$fa_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]-1/3,-abs(anagenetic_fuzzies$la_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+1/3,col=makeTransparent(otu_colors[sp],alpha=100),border=makeTransparent(otu_colors[sp],alpha=100),lwd=0);
#			if (abs(anagenetic_fuzzies$fa_ub[fuzzy_dude])>abs(anagenetic_fuzzies$la_lb[fuzzy_dude]))
#				rect(-abs(anagenetic_fuzzies$fa_ub[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]-1/3,-abs(anagenetic_fuzzies$la_lb[fuzzy_dude]),anagenetic_fuzzies$phylo_axis[fuzzy_dude]+1/3,col=otu_colors[sp],lwd=0);
#		if ()
			}
		}
	}
observed_nodes_all <- rep(0,nNodes);
an <- 0;
# draw out ancestors specially
while (an < length(ancestral_spc))	{
	an <- an+1;
	observed_nodes_all[which(mat_tree==ancestral_spc[an],arr.ind = T)[1]] <- ancestral_spc[an];
	}
# draw branches linking taxa to nodes and nodes to deeper nodes
for (nn in 1:length(phylo_axis))	{
#	nn <- nn+1;
	if (nn<=notu)	{
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
# link nodes & branches
for (nd in 1:nNodes)	{
	htu <- nd+notu;
	f1 <- mat_tree[nd,mat_tree[nd,]>0];
	segments(branch_ranges$finish[htu],min(phylo_axis[f1]),branch_ranges$finish[htu],max(phylo_axis[f1]));
	}
# draw reconstructed FAs
for (sp in 1:notu)	{
	points(-abs(branch_ranges$finish[sp]),phylo_axis[sp],pch=8,cex=1.0,lwd=1.5);
	points(-abs(branch_ranges$finish[sp]),phylo_axis[sp],pch=8,cex=1.0,lwd=2/3,col=otu_colors[sp]);
	}
top <- mxy+1;
taxon_font_size <- 0.6;
if (print_names)	{
	sp <- 0;
#	for (sp in 1:notu)	{
	while (sp < notu)	{
		sp <- sp+1;
		if (!sp %in% anagenetic_fuzzies$anc_no && !sp %in% anagenetic_fuzzies$dsc_no)	{
			xx <- max(branch_ranges$finish[sp],-abs(fuzzy_ranges$la_ub[sp]));
			text(xx,phylo_axis[sp],otu_names[sp],pos=4,cex=taxon_font_size,family=franky);
			} else	if (sp %in% anagenetic_fuzzies$anc_no)	{
			fuzzy_dude <- match(sp,anagenetic_fuzzies$anc_no);
			xx <- max(branch_ranges$finish[sp],-anagenetic_fuzzies$la_ub[fuzzy_dude]);
#		lineage_name <- paste(anagenetic_fuzzies$ancestor[fuzzy_dude],"->",anagenetic_fuzzies$descendant[fuzzy_dude],sep="")
#			xx <- max(branch_ranges$finish[sp],-abs(fuzzy_ranges$la_ub[sp]));
			text(xx,phylo_axis[sp],anagenetic_fuzzies$taxon[fuzzy_dude],pos=4,cex=taxon_font_size,family=franky);
			} else	if (sp %in% anagenetic_fuzzies$dsc_no)	{
#			fuzzy_dude <- match(sp,anagenetic_fuzzies$dsc_no);
#		lineage_name <- paste(anagenetic_fuzzies$ancestor[fuzzy_dude],"->",anagenetic_fuzzies$descendant[fuzzy_dude],sep="")
#			xx <- max(branch_ranges$finish[sp],-abs(fuzzy_ranges$la_ub[sp]));
#			text(xx,phylo_axis[sp],anagenetic_fuzzies$taxon[fuzzy_dude],pos=4,cex=taxon_font_size,family=franky);
			}
		}
	}

# print separate key; ####
# taxon groups:
group_color_scheme <- unique(otu_information[,c(rank_to_plot,"color","taxon_type")]);
group_color_scheme <- group_color_scheme[group_color_scheme$taxon_type %in% "ingroup",];
#if (!"outgroup" %in% tolower(group_color_scheme[,rank_to_plot]))
#	group_color_scheme <- rbind(group_color_scheme,c("outgroup","gray75"));
group_colors <- group_color_scheme$color;
names(group_colors) <- group_color_scheme[,rank_to_plot];
x_to_y <- ((end-onset)/xsize)/(mxy/ysize);
y_to_x <- (mxy/ysize)/((end-onset)/xsize);

#rect(onset+(x_to_y*0),10,onset+(x_to_y*2),8);
#rect(onset+(x_to_y*2),10,onset+(x_to_y*4),8);


group_colors <- group_colors[match(unique(group_colors),group_colors)];
group_colors <- group_colors[order(names(group_colors))]
group_colors <- group_colors[!group_colors %in% "pink3"]
if (tree_to_print %in% c("1L","1U"))	group_colors <- group_colors[match(c("Asteroblastida","Glyptosphaeritida","Sphaeronitida","Platycystitida","Rhombifera","Crinoidea","Coronata","Parablastoidea"),names(group_colors),)];
specify_basic_plot(mxx=2*length(group_colors),mnx=0,mxy=2*length(group_colors),mny=0,xsize=2,ysize=2);
for (i in length(group_colors):1)	{
	j <- 1+length(group_colors)-i;
	rect(0.0,i+0.9,0.8,i+0.1,col="white",border=group_colors[j]);
	rect(0.0,i+0.9,0.4,i+0.1,col=makeTransparent(group_colors[j],alpha=100),lwd=0);
	rect(0.4,i+0.9,0.8,i+0.1,col=group_colors[j],lwd=0);
	text(0.5,i+0.5,paste(":",names(group_colors)[j]),pos=4,cex=1.25*otu_font_size,family=franky);
	}

specify_basic_plot(mxx=6,mnx=0,mxy=12,mny=0,xsize=2,ysize=2);
for (i in 5:1)	{
	segments(0,i+0.5,1,i+0.5,lwd=5*i/5);
	wibble <- paste(": posterior prob. =",round(i/5,1));
	if (i==5)	wibble <- paste(wibble,".0",sep="");
	text(0.8,i+0.5,wibble,pos=4,cex=1.25*otu_font_size,family=franky);
	}
