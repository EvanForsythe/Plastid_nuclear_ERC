###Script for running full ERC pipeline.
#The pipeline takes nuclear gene trees and filters/processes them to be suitable for ERC analysis
#ERC is then performed comparing branch lengths from these trees against branch lengths for plastid trees (also input)
#The final output is a csv file (for each plastid partition) with correlation statistics for each nuclear gene tree

#Name working directory (don't forget to change to match your local environment)
working_dir<-"/Plastid_nuclear_ERC/"

#Set working directory
setwd(working_dir)

##Load packages
package_list<-c("plyr", "seqinr", "ggplot2", "gplots", "RColorBrewer", "ape", "insect", "phytools", "ggplot2", "grid", "treeio", "phangorn")

#Loop to check if package is installed and libraried
for(p in 1:length(package_list)){
  if (!require(package_list[p], character.only = TRUE)) {
    install.packages(package_list[p], dependencies = TRUE)
    library(package_list[p], character.only=TRUE)
  }
}

#Install biomanager in order to install treeio
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

#Install and library treeio
#BiocManager::install("treeio")
library(treeio)

################################################
#######  Process the full gene trees  ##########
#######  and extract mininum subtrees ##########
################################################
#Make lists of species names
All_eudicot<-c(
  "O_biennis","E_grandis","C_sativus","P_persica","A_aulacocarpa",
  "G_raimondii","P_trichocarpa","V_vinifera","G_maderense","A_thaliana",
  "L_siphilitica","H_annuus","P_maritima","S_lycopersicum", "S_noctiflora"
)

All_monocots<-c(
  "M_acuminata","O_sativa","S_polyrhiza"
)

All_outgroup<-c("A_trichopoda", "L_chinense")

#Name the directory that the trees that will be input
tree_dir<-paste(working_dir,"Trees/", sep = "")

#Get list of all input trees in the directory
In_tree_files<-list.files(path = tree_dir, pattern = "RAxML_bipartitions*")

# #It may necissary to exclude rearranged trees (from previous runs)
# if(length(In_tree_files[-grep("rearrange", In_tree_files)])>0){
# In_tree_files<-In_tree_files[-grep("rearrange", In_tree_files)]
# }

#Make multiphylo to store all the subtrees that pass filters (keeper subtrees)
Keeper_subtrees<-list()
class(Keeper_subtrees)<-"multiPhylo"

#Make a matrix to store descriptive statistics for the trees
tree_stats_mat<-matrix(, ncol = 3, nrow = length(In_tree_files))

### Loop through all trees ###
for(t in 1:length(In_tree_files)){

#Rearrange via GT/SP reconciliation in Notung by calling Notung from the command line
#make the command
notung_rearrange_cmd<-paste("java -jar /Applications/Notung-2.9/Notung-2.9.jar ", tree_dir, In_tree_files[t], " -s SP_tree.newick --rearrange --threshold 80 --treeoutput nhx --nolosses --speciestag prefix --edgeweights name --outputdir ", tree_dir, sep = "")
#Run the command
system(notung_rearrange_cmd)

#Root in Notung
#make the command
notung_root_cmd<-paste("java -jar /Applications/Notung-2.9/Notung-2.9.jar ", tree_dir, In_tree_files[t], ".rearrange.0 ", " -s SP_tree.newick --root --treeoutput nhx --nolosses --speciestag prefix --edgeweights name --outputdir ", tree_dir, sep = "")
#Run the command
system(notung_root_cmd)

#Read in the GT/ST reconciled tree
tree<-read.tree(file = paste(tree_dir, In_tree_files[t], ".rearrange.0.rooting.0", sep = ""))

#Read in the same tree but output in NHX format
nhx_tree<-read.nhx(file = paste(tree_dir, In_tree_files[t], ".rearrange.0.rooting.0", sep = ""))

#Set all branch lengths to 1 (we're just using topology at this point)
tree$edge.length<-sample(1, length(tree$edge.length), replace = TRUE)

#Get dup nodes from the nhx
dup_nodes<-grep("Y", nhx_tree@data$D)

#Make multiphylo to store clades
full_clades<-list()
class(full_clades)<-"multiPhylo"

#Loop through each of the duplication nodes
for(n in 1:length(dup_nodes)){
  #Get the two decendent nodes from the dup nodes
  des_nodes<-Children(tree, dup_nodes[n])
  
  #Loop through each of the two decendent nodes
  for(d in 1:length(des_nodes)){
    #Ask if extract clade throws an error, proceed if not
    if(length(grep("Error", try(extract.clade(tree, des_nodes[d]), silent = TRUE))) == 0){
  
      #Extract just the subtree from the decendent node  
  temp_clade<-extract.clade(tree, des_nodes[d]) 

  #Make lists to to store tip counts
  e_tips<-list()
  m_tips<-list()
  o_tips<-list()
  
  #Make list of eudicot tips
  for(e in 1:length(All_eudicot)){
    e_tips<-c(e_tips, temp_clade$tip.label[grep(All_eudicot[e], temp_clade$tip.label)])
  }
  #Make list of monocot tips
  for(m in 1:length(All_monocots)){
    m_tips<-c(m_tips, temp_clade$tip.label[grep(All_monocots[m], temp_clade$tip.label)])
  }
  #Get outgroup tips
  for(o in 1:length(All_outgroup)){
    o_tips<-c(o_tips, temp_clade$tip.label[grep(All_outgroup[o], temp_clade$tip.label)])
  }
  
  #Check if subtree contains minimum taxon sampling and monophyly requirements
  if(length(e_tips) > 0 & length(m_tips) > 0 
     & length(o_tips) > 0 & length(unique(substr(c(e_tips, m_tips), 1, 5))) > 9
     & is.monophyletic(temp_clade, unlist(c(e_tips, m_tips)))
     ){
      
    #Add keeper to keeper list
    full_clades[[length(full_clades)+1]]<-temp_clade
    
      } #End taxon composition if statement
    } #End if statement checking for error
  } #End descendent nodes loop
} #End dup nodes loop

#If no subtrees were retained above, check out the full tree to see if it should be retained
if(length(full_clades)==0){
#Make places to store tip counts
e_tips2<-list()
m_tips2<-list()
o_tips2<-list()

#Make list of eudicot tips
for(e in 1:length(All_eudicot)){
  e_tips2<-c(e_tips2, tree$tip.label[grep(All_eudicot[e], tree$tip.label)])
}
#Make list of monocot tips
for(m in 1:length(All_monocots)){
  m_tips2<-c(m_tips2, tree$tip.label[grep(All_monocots[m], tree$tip.label)])
}
#Get outgroup tips
for(o in 1:length(All_outgroup)){
  o_tips2<-c(o_tips2, tree$tip.label[grep(All_outgroup[o], tree$tip.label)])
}

#Check if tree contains minimum taxon sampling and monophyly requirements
if(length(e_tips2) > 0 & length(m_tips2) > 0 
   & length(o_tips2) > 0 & length(unique(substr(c(e_tips2, m_tips2), 1, 5))) > 9
   & is.monophyletic(tree, unlist(c(e_tips2, m_tips2)))
   ){
  
  #Add full trees to list
  full_clades[[length(full_clades)+1]]<-tree
  
}#End taxon composition loop
}#End full tree if statement

#Check if there are any full trees
#If there are, we need to go through them to find the clades that don't overlap with eachother
if(length(full_clades)>0){

#Make multiphylo object to store clades
full_nonoverlap_clades<-list()
class(full_nonoverlap_clades)<-"multiPhylo"

#Loop through the clades
Keeper_taxa<-list()
while(length(full_clades) > 0){
  
  #Get taxon counts for trees in full_clades
  #Make empty list for taxon counts
  tax_counts<-list()
  
  #Loop through the full clades to check for overlap
  for(k in 1:length(full_clades)){
    tax_counts[k]<-length(full_clades[[k]]$tip.label)
  }
  
  #Get the smallest clade
  small_clade<-full_clades[[which(tax_counts == min(unlist(tax_counts)))[1]]]
  
  #Get small clade taxa
  small_clade_taxa<-small_clade$tip.label
  
  #Start an emplty list to store overlap taxa
  overlap_hits<-list()
 
  #Loop through small clade taxa
  for(s in 1:length(small_clade_taxa)){
    overlap_hits[s]<-length(grep(small_clade_taxa[s], Keeper_taxa))
  }
  
  #Ask if there's overlap
  if(!any(unlist(overlap_hits) > 0)){
    #Add to full+non-overlapping list if it's not overlapping
    full_nonoverlap_clades[[length(full_nonoverlap_clades)+1]]<-small_clade
    Keeper_taxa<-unlist(c(Keeper_taxa, small_clade$tip.label))
  }
  
  #Make a new "full_clades" with all clades except for the previous small_clade
  remaining_clades<-seq(1, length(full_clades))[-which(tax_counts == min(unlist(tax_counts)))[1]]
  
  #Go through each of the remaining clades to populate a new multiphylo  
  #Make multiphylo to store new list of clades (omiting small_clade)
  full_clades_temp<-list()
  class(full_clades_temp)<-"multiPhylo"
  
  #Check to see if there are any remaining clades
  if(length(remaining_clades)>0){
  #Loop through all the remaining trees
  for(r in 1:length(remaining_clades)){
      full_clades_temp[[r]]<-full_clades[[remaining_clades[r]]]
    }#End for statement that builds full_clades_temp
  }#End remaining clades if statement
  
  #Move full_clades_temp (smallest clade removed) to full_clades
  full_clades<-full_clades_temp

} #End while loop

#Add the current keepers to the full keeper collection
for(s in 1:length(full_nonoverlap_clades)){
Keeper_subtrees[[length(Keeper_subtrees)+1]]<-full_nonoverlap_clades[[s]]
}

#Fill in the stats
tree_stats_mat[t,]<-c(In_tree_files[t], length(tree$tip.label), length(full_nonoverlap_clades))
}else{
  tree_stats_mat[t,]<-c(In_tree_files[t], length(tree$tip.label), 0)
}#End if/else to see if there are any full clades
}#End all trees for loop

## Write results of subtree process
write.csv(tree_stats_mat, file = "tree_stats_mat.csv")

#Write each of the keeper subtrees as a newick file
for(n in 1:length(Keeper_subtrees)){
write.tree(Keeper_subtrees[[n]], file = paste(working_dir, "Subtrees/Subtree_", sprintf("%04d", n), ".txt", sep = ""))
}

### Make plots of descriptive stats about subtree pruning
#Convert to dataframe
tree_stats_df<-data.frame(File_name=paste(tree_stats_mat[,1]), 
           Taxon_counts=as.numeric(paste(tree_stats_mat[,2])),
           Retained_subtrees=as.numeric(paste(tree_stats_mat[,3]))
           )

#Convert NAs to 0s for plotting
tree_stats_df$Taxon_counts[which(is.na(tree_stats_df$Taxon_counts)==TRUE)]<-0

#Make a quick plot of retained subtrees vs taxa in full tree
plot(y=tree_stats_df$Taxon_counts, x=jitter(tree_stats_df$Retained_subtrees),
     xlab = "Retained subtrees", ylab = "Taxa in full tree")

#Make histogram of number taxa in subtrees
taxon_counts<-list()
for(s in 1:length(Keeper_subtrees)){
  taxon_counts[s]<-length(Keeper_subtrees[[s]]$tip.label)
}

#Plot the histogram
hist(unlist(taxon_counts), breaks = 100, xlab = "Total taxa/subtree")
abline(v = 10)

#Make histogram of species represented in subtrees
species_counts<-list()
#p<-1
for(p in 1:length(Keeper_subtrees)){
  species_counts[p]<-length(unique(substr(Keeper_subtrees[[p]]$tip.label, 1, 5)))
}

#Plot the histogram
hist(unlist(species_counts), xlab = "Species represented/subtree", breaks = 10)

#If "Keeper_subtrees" object doesn't exist, create it by reading in subtrees 
if(!exists("Keeper_subtrees")){

subtree_list<-list.files(path = paste(working_dir, "Subtrees/", sep = ""), pattern = "Subtree_")
Keeper_subtrees<-list()
class(Keeper_subtrees)<-"multiPhylo"
for(k in 1:length(subtree_list)){
  Keeper_subtrees[[k]]<-read.tree(paste(working_dir, "Subtrees/", subtree_list[k], sep = ""))
}
}#End if statement

#################################
#######  Subtree parsing  #######
#################################
#Now that we've identified the subtrees to work with, we need to pull out
#the sequences from the FASTA files so we can realign and reinfer trees

#Extract sequence from subtrees
#Define the directory to work in
subtree_parse_dir<-paste(working_dir,"Subtree_seq_parsing/", sep = "")

#Make blank matrix
Keeper_trees_mat<-matrix(, ncol = 110, nrow = length(Keeper_subtrees))

#Fill in matrix
for(s in 1:length(Keeper_subtrees)){
  Keeper_trees_mat[s,1]<-paste("Subtree_", sprintf("%04d", s), sep = "")
  Keeper_trees_mat[s,2:(length(Keeper_subtrees[[s]]$tip.label)+1)]<-c(Keeper_subtrees[[s]]$tip.label)
  write.tree(Keeper_subtrees[[s]], file = paste(subtree_parse_dir, "Raw_subtrees/","Subtree_", sprintf("%04d", s), sep = ""))
}

#Write the csv file that contains the seq IDs in each subtrees
write.csv(Keeper_trees_mat, file = paste(subtree_parse_dir, "Min_subtrees_forERC.csv", sep = ""), row.names = FALSE, col.names = FALSE)

#Extract subtree fastas from large fasta by calling a python script
extract_cmd<-paste("python3.6 ", subtree_parse_dir, "makeSeqFiles.py ", subtree_parse_dir, "Min_subtrees_forERC.csv ", subtree_parse_dir, "ALL_PROTS.fa ", sep = "")

#Run command
system(extract_cmd)

#Run these commands to move files to the directory the belong to
#This is just some quick housekeeping
#The mv commands need to be done in a few steps because there's so many files it throws an error.
system("mv Subtree_1*_seqs.fasta Subtree_seq_parsing/")
system("mv Subtree_2*_seqs.fasta Subtree_seq_parsing/")
system("mv Subtree_*_seqs.fasta Subtree_seq_parsing/")

# Perform alignment of the new subtrees
subtree_seqs_list<-list.files(subtree_parse_dir, pattern = "Subtree")

#name the directory for the alignmnts
New_aln_dir<-paste(working_dir,"Aligned_subtrees/", sep = "")

#Loop through all the subtrees and run mafft command
for(s in 1:length(subtree_seqs_list)){
  mafft_cmd<-paste("mafft-linsi --thread 6 ", subtree_parse_dir, subtree_seqs_list[s], " > ", New_aln_dir, "aln_",subtree_seqs_list[s], sep = "")
  #Run the command
  system(mafft_cmd)
}

###############################################
########  Clean alignment with gblocks ########
###############################################
#DSownloaded from: http://molevol.cmima.csic.es/castresana/Gblocks.html

#Get the names of all the files in the alignent folder
alns<-list.files(path = New_aln_dir, pattern = "aln_")

#Trim this list to exclude already gblocked alignments (in case this block of code is repeated for any reason)
if(length(grep("-gb", alns))>0){
alns<-alns[-(grep("-gb", alns))]
}

#Loop through all the trees
for(a in 1:length(alns)){
  #Create command to count the number of seqs in the aln
  grep_cmd<-paste("grep '>' ", New_aln_dir, alns[a], " | wc -l", sep = "")
  num_seqs<-gsub("[[:blank:]]", "", system(grep_cmd, intern = TRUE))
  
  #Convert this to the desired number for the parameter
  num_seq_param<-floor((as.numeric(num_seqs)/2)+1)
  
  #Create command for running gblocks
  gblocks_cmd<-paste("Gblocks ", New_aln_dir, alns[a], " -b5=h -b4=5 -b2=", num_seq_param, sep = "")
  #Run gblocks
  system(gblocks_cmd)

  }#End aln for-loop

### Next, I build trees from the alignments on a super computer (below are examples of the executables I used)
#command for a batch script
#"raxmlHPC-PTHREADS-SSE3 -s $file -n $file -m PROTGAMMALGF -p 12345 -x 12345 -# 100 -f a -T 24"
# If doing locally for some reason: use "raxmlHPC-PTHREADS-AVX"

### Inferred branch lengths on the reinferred trees with the follwing command:
#"raxmlHPC-PTHREADS-SSE3 -s aln_${file}_seqs.fasta -n ${file} -t ${file} -m PROTGAMMALGF -p 12345 -T 24 -f e"

#Work with the branch length optimized trees to root by the root obtained from the subtree extraction above
### Get dir where trees live
raw_subtrees_dir<-paste(working_dir, "Mono_required_80bs/Subtree_seq_parsing/Raw_subtrees/", sep = "")
BL_subtrees_dir<-paste(working_dir, "Mono_required_80bs/BL_trees/", sep = "")

raw_subtrees_list<-list.files(path = raw_subtrees_dir, pattern = "Subtree_")

#List to store rooted trees
nuc_BL_trees<-list()
class(nuc_BL_trees)<-"multiPhylo"

#r<-1
#Loop through all raw subtrees
for(r in 1:length(raw_subtrees_list)){

#read in the tree
raw_tree<-read.tree(file = paste(raw_subtrees_dir, raw_subtrees_list[r], sep = ""))
BL_tree<-read.tree(file = paste(BL_subtrees_dir, "RAxML_result.", raw_subtrees_list[r], sep = ""))

subtrees<-treeSlice(raw_tree, 0.01, trivial=FALSE, prompt=FALSE)

if(length(subtrees)==1){
  outtaxa<-raw_tree$tip.label[!raw_tree$tip.label %in% subtrees[[1]]$tip.label]
}else if(length(subtrees)==2){

t_one_count<-subtrees[[1]]$tip.label
t_two_count<-subtrees[[2]]$tip.label

if(length(t_one_count)<length(t_two_count)){
  outtaxa<-t_one_count
}else if(length(t_one_count)>length(t_two_count)){
  outtaxa<-t_two_count
}
} #End subtrees if statement.

#Root the tree
raw_tree_root<-root(BL_tree, outgroup = outtaxa, resolve.root = TRUE)

#Add the rooted tree to the multiphylo
nuc_BL_trees[[r]]<-raw_tree_root

}#End all trees loop

####################################
#### WORK WITH REINFERRED TREES ####
####################################
reinfer_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC/ERC_PIPELINE/Trim_Retree/Reconciled_subtrees/"

#Get file names
reconciled_files<-list.files(path = reinfer_dir, pattern = "rearrange.0.rooting.0")

for(c in 1:length(reconciled_files)){
  tree<-read.tree(file = paste(reinfer_dir, reconciled_files[c], sep = ""))
  tree$node.label<-NULL
  #write.tree(tree, file = substr(reconciled_files[c], 24, 35))
}

#####################################################
########    Filter trees before ERC    ##############
#####################################################

#Get directory where trees and alignments live
gb_aln_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC/ERC_PIPELINE/Trim_Retree/Trimmed_alns/"
bl_tree_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC/ERC_PIPELINE/Trim_Retree/BL_trees/"

#Make a range of numbers cooresponding to the subtree numbers. We have 8150 subtrees
subtree_nums<-seq(1, 8150, by = 1)

#make a blank list to store the names of alignments that are too short (alignment length shitlist)
aln_ln_shitlist<-list()

#Loop through all alignments to look at length
for(a in 1:length(subtree_nums)){
  aln_temp<-read.alignment(file = paste(gb_aln_dir, "aln_Subtree_", sprintf("%04d", subtree_nums[a]), "_seqs.fasta-gb", sep = ""), format = "fasta")
  if(nchar(aln_temp$seq[1])<50){
    aln_ln_shitlist[length(aln_ln_shitlist)+1]<-sprintf("%04d", subtree_nums[a])
  }
}

#Make list for names of trees with long outlier branches that we want to throw out (outlier shitlist)
tree_outlier_shitlist<-list()
tree_outlier_ratio<-list()
full_prune_shitlist<-list()

#Look through the longest and second longest branches
for(t in 1:length(subtree_nums)){
  if(file.exists(file = paste(bl_tree_dir, "RAxML_result.Subtree_", sprintf("%04d", subtree_nums[t]), sep = ""))){
  tree_temp<-read.tree(file = paste(bl_tree_dir, "RAxML_result.Subtree_", sprintf("%04d", subtree_nums[t]), sep = ""))
  tree_temp<-unroot(tree_temp)
  longest<-tree_temp$edge.length[order(-tree_temp$edge.length)][1]
  second_longest<-tree_temp$edge.length[order(-tree_temp$edge.length)][2]
  tree_outlier_ratio[t]<-longest/second_longest
  if(longest/second_longest>10){
    tree_outlier_shitlist[length(tree_outlier_shitlist)+1]<-sprintf("%04d", subtree_nums[t])
  }#End local if statement
  }else{
    full_prune_shitlist[length(full_prune_shitlist)+1]<-sprintf("%04d", subtree_nums[t])
  }
}

#Make a quick venn diagram of the three shitlists
venn_shitlists<-list(unlist(full_prune_shitlist), unlist(tree_outlier_shitlist), unlist(aln_ln_shitlist))
venn::venn(venn_shitlists, snames = c("seq_lost_gblocks", "branch_length", "aln_length"), zcolor = brewer.pal(n = 3, name = "Accent"))

#Get the final full shitlist
FULL_shitlist<-unique(c(unlist(full_prune_shitlist), unlist(tree_outlier_shitlist), unlist(aln_ln_shitlist)))

#length(FULL_shitlist)/8150
#220 genes were lost from all three filters combined
#This is 2.69% of the subtrees

#Get the final list of trees that will be input into ERC
final_keeper_subtrees<-subtree_nums[-as.numeric(FULL_shitlist)]

#Read in each tree and root by the outgroup from notung
noBL_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC/ERC_PIPELINE/Trim_Retree/Reconciled_subtrees/"
BL_dir<-"/Users/esforsythe/Documents/Work/Bioinformatics/ERC/ERC_PIPELINE/Trim_Retree/BL_trees/"

#List to store rooted trees
nuc_BL_trees<-list()
class(nuc_BL_trees)<-"multiPhylo"

# Do a quick test of monophyly 
Monophyly_test<-list()

for(k in 1:length(final_keeper_subtrees)){
  raw_tree<-read.tree(file = paste(noBL_dir, "Subtree_", sprintf("%04d", final_keeper_subtrees[k]), sep = ""))
  #Set all branch lengths on raw tree to 1 so that treeSlice will work
  raw_tree$edge.length[1:length(raw_tree$edge.length)]<-1
  bl_tree<-read.tree(file = paste(BL_dir, "RAxML_result.Subtree_", sprintf("%04d", final_keeper_subtrees[k]), sep = ""))

  #Get the outgroup taxa for raw tree
  Amborella_seqs<-raw_tree$tip.label[grep("trichopoda" ,raw_tree$tip.label)]
  Lirio_seqs<-raw_tree$tip.label[grep("chinense" ,raw_tree$tip.label)]
  full_outgroup<-raw_tree$tip.label[c(grep("trichopoda" ,raw_tree$tip.label), grep("chinense" ,raw_tree$tip.label))]
  full_ingroup<-raw_tree$tip.label[-c(grep("trichopoda" ,raw_tree$tip.label), grep("chinense" ,raw_tree$tip.label))]
  
  #Slice tree at root
  subtrees<-treeSlice(raw_tree, 0.01, trivial=FALSE, prompt=FALSE)
    
    if(length(subtrees)==1){
      outtaxa<-raw_tree$tip.label[!raw_tree$tip.label %in% subtrees[[1]]$tip.label]
    }else if(length(subtrees)==2){
      
      t_one_count<-subtrees[[1]]$tip.label
      t_two_count<-subtrees[[2]]$tip.label
      
      if(length(t_one_count)<length(t_two_count)){
        outtaxa<-t_one_count
        intaxa<-t_two_count
      }else if(length(t_one_count)>=length(t_two_count)){
        outtaxa<-t_two_count
        intaxa<-t_one_count
      }
    } #End subtrees if statement.
    
    ##### Go through the trees to and prune outgroup branches if they are nested within the ingroup and reroot by the remaining outgroup(s)
    #Ask if ingroup is monophyletic
    if(!is.monophyletic(raw_tree, full_ingroup)){ #If it's not monophyletic, do the following
      Monophyly_test[k]<-FALSE
      if(length(intersect(outtaxa, full_outgroup))==length(outtaxa)){ #True if all rooting taxa are outtaxa
        bl_tree_root<-root(bl_tree, outgroup = outtaxa, resolve.root = TRUE)
        bl_tree_root<-drop.tip(bl_tree_root, tip = 
                                full_outgroup[which(is.element(full_outgroup, outtaxa)==FALSE)]
                                ) #End drop tip command
          }else if(!length(intersect(outtaxa, full_outgroup))==length(outtaxa)){ #True if NOT all rooting taxa are outtaxa
            #choose a random Amoborella (or Liriodenedron if Amberella isnt there) and prune (drop) all others
            if(length(Amborella_seqs)>0){ #if there is amborella
            bl_tree_root<-root(bl_tree, outgroup = Amborella_seqs[1], resolve.root = TRUE)
            bl_tree_root<-drop.tip(bl_tree_root, tip = c(Amborella_seqs[2:length(Amborella_seqs)], Lirio_seqs))
          }else{ #else if there's no amborella
            bl_tree_root<-root(bl_tree, outgroup = Lirio_seqs[1], resolve.root = TRUE)
            bl_tree_root<-drop.tip(bl_tree_root, tip = Lirio_seqs[2:length(Lirio_seqs)])
          }# end no amborella else
          }# end all outgroup not outtaxa else
    }else if(is.monophyletic(raw_tree, full_ingroup)){ #If ingroup is monophyletic
      Monophyly_test[k]<-TRUE
      bl_tree_root<-root(bl_tree, outgroup = outtaxa, resolve.root = TRUE)
    } #End is monophyletic if
    
    #Add the rooted tree to the multiphylo
    nuc_BL_trees[[length(nuc_BL_trees)+1]]<-bl_tree_root
  } #End keeper trees loops

### Write a mulitphlyo file of nuclear BL trees
write.tree(nuc_BL_trees, file = "nuclear_rooted_BL_trees.txt")

#####################################################
##########  GET BRANCH LENGTHS FOR ERC  #############
#####################################################

#Get Plastid trees
plastid_trees<-read.tree(file = "Plastid_BL_trees_cat")

#Get the nuclear trees if the object doesn't exist
nuc_BL_trees<-read.tree(file="nuclear_rooted_BL_trees.txt")

rooted_plast<-list()
class(rooted_plast)<-"multiPhylo"

#Root the plastid trees
for(t in 1:length(plastid_trees)){
  rooted_plast[[t]]<-root(plastid_trees[[t]], outgroup = which(plastid_trees[[t]]$tip.label=="Amborella_trichopoda"), resolve.root = TRUE)
}

#Make list of species in the tree (Note: I excluded "chinense" and "trichopoda" from this list because)
sp_list<-c("raimondii", "sativus", "noctiflora", "vinifera", "lycopersicum", "sativa",
           "acuminata", "trichocarpa", "grandis", "biennis",
           "thaliana", "siphilitica", "persica", "maritima", "maderense", "aulacocarpa", 
           "polyrhiza", "annuus")

###Combine plastid and nuclear trees in order to get branch lengths from them all
nuc_and_cp_trees<-c(rooted_plast, nuc_BL_trees)

#Make a list of plastid gene names
plast_names<-c("RNApol", "ACCD", "CLPP1", "MATK", "PHOTOSYNTH", "RIBO", "YCF")

#Make matrix to store trees results
dist_out_g_trees<-matrix(, ncol = length(sp_list), nrow = length(nuc_and_cp_trees))

#Make blank lists
Athal_names<-list()
Osat_names<-list()
Subtree_file<-list()
Monophyly_test_new<-list()

#Loop throgh all the trees
for(g in 1: length(nuc_and_cp_trees)){

  #Get one tree from the list
  tree<-nuc_and_cp_trees[[g]]
  
  #Get rosids taxa in the tree
  All_rosids<-c(
    grep("biennis", tree$tip.label),
    grep("grandis", tree$tip.label),
    grep("sativus", tree$tip.label),
    grep("persica", tree$tip.label),
    grep("aulacocarpa", tree$tip.label),
    grep("raimondii", tree$tip.label),
    grep("trichocarpa", tree$tip.label),
    grep("vinifera", tree$tip.label),
    grep("maderense", tree$tip.label),
    grep("thaliana", tree$tip.label)
  )
  
  #Get asterid taxa in the tree
  All_asterids<-c(
    grep("siphilitica", tree$tip.label),
    grep("annuus", tree$tip.label),
    grep("maritima", tree$tip.label),
    grep("lycopersicum", tree$tip.label),
    grep("noctiflora", tree$tip.label)
  )
  
  #Get monocot taxa in the tree
  All_monocots<-c(
    grep("acuminata", tree$tip.label),
    grep("sativa", tree$tip.label),
    grep("polyrhiza", tree$tip.label)
  )
  
  #Get the node uniting ingroup
  ingroup_nodes<-mrca.phylo(tree, c(All_rosids, All_asterids, All_monocots))
  
  #Check to make sure ingroup is monophyletic
  Monophyly_test_new[g]<-is.monophyletic(tree, c(All_rosids, All_asterids, All_monocots))
  
  #Get the distance of each taxon to the ingroup node
  ingroup_distances<-dist.nodes(tree)[ingroup_nodes,]
  
  #Make blank list to store instances 
  sp_dist_out<-list()
  
  #Extract distance for each species (take the average if there's paralogs)
  for(n in 1: length(sp_list)){
    sp_dist_out[n]<-mean(ingroup_distances[grep(sp_list[n], tree$tip.label)])
  }
  
  #Output distances and other info
  dist_out_g_trees[g,]<-unlist(sp_dist_out) 
  Athal_names[g]<-paste(tree$tip.label[grep("thaliana", tree$tip.label)], collapse = "")
  Osat_names[g]<-paste(tree$tip.label[grep("sativa", tree$tip.label)], collapse = "")
  
  #Do this only for the nuclear trees (not the plastid trees at the top of the list)
  if(g>length(rooted_plast)){
    #Subtree_file[g]<-raw_subtrees_list[g-length(rooted_plast)]
    Subtree_file[g]<-final_keeper_subtrees[g-length(rooted_plast)]
  }
}#End loop of all trees (nuc and plast)

#####################################################
##########  CLEAN BRANCH-LENGTH RESULTS #############
#####################################################

#convert to dataframe and rename 
dist_out_g_trees_df<-as.data.frame(dist_out_g_trees)
names(dist_out_g_trees_df)<-sp_list

#Make a blank matrix to store root-to-tip branclengths
root2tip_stats<-matrix(, nrow= ncol(dist_out_g_trees_df), ncol=3)

### Take a look at the distributions of branch lengths (root2tips)
par(mfrow=c(6,3))
for(r in 1:ncol(dist_out_g_trees_df)){
  plot(density(na.omit(dist_out_g_trees_df[,r])), xlim = c(0,5), main = names(dist_out_g_trees_df)[r], xlab = NULL)
  abline(v=mean(na.omit(dist_out_g_trees_df[,r])))
  abline(v=median(na.omit(dist_out_g_trees_df[,r])), col = "red", lty =2)

  root2tip_stats[r,]<-c(names(dist_out_g_trees_df)[r], mean(na.omit(dist_out_g_trees_df[,r])), median(na.omit(dist_out_g_trees_df[,r])))

}
par(mfrow=c(1,1))

#Normalize all values by the genome-wide MEAN
dist_out_norm_df<-data.frame(
  G_raimondii=dist_out_g_trees_df$raimondii/mean(na.omit(dist_out_g_trees_df$raimondii)),
  C_sativus=dist_out_g_trees_df$sativus/mean(na.omit(dist_out_g_trees_df$sativus)),
  S_noctiflora=dist_out_g_trees_df$noctiflora/mean(na.omit(dist_out_g_trees_df$noctiflora)),
  V_vinifera=dist_out_g_trees_df$vinifera/mean(na.omit(dist_out_g_trees_df$vinifera)),
  S_lycopersicum=dist_out_g_trees_df$lycopersicum/mean(na.omit(dist_out_g_trees_df$lycopersicum)),
  O_sativa=dist_out_g_trees_df$sativa/mean(na.omit(dist_out_g_trees_df$sativa)),
  M_acuminata=dist_out_g_trees_df$acuminata/mean(na.omit(dist_out_g_trees_df$acuminata)),
  P_trichocarpa=dist_out_g_trees_df$trichocarpa/mean(na.omit(dist_out_g_trees_df$trichocarpa)),
  E_grandis=dist_out_g_trees_df$grandis/mean(na.omit(dist_out_g_trees_df$grandis)),
  O_biennis=dist_out_g_trees_df$biennis/mean(na.omit(dist_out_g_trees_df$biennis)),
  A_thaliana=dist_out_g_trees_df$thaliana/mean(na.omit(dist_out_g_trees_df$thaliana)),
  L_siphilitica=dist_out_g_trees_df$siphilitica/mean(na.omit(dist_out_g_trees_df$siphilitica)),
  P_persica=dist_out_g_trees_df$persica/mean(na.omit(dist_out_g_trees_df$persica)),
  P_maritima=dist_out_g_trees_df$maritima/mean(na.omit(dist_out_g_trees_df$maritima)),
  G_maderense=dist_out_g_trees_df$maderense/mean(na.omit(dist_out_g_trees_df$maderense)),
  A_aulacocarpa=dist_out_g_trees_df$aulacocarpa/mean(na.omit(dist_out_g_trees_df$aulacocarpa)),
  S_polyrhiza=dist_out_g_trees_df$polyrhiza/mean(na.omit(dist_out_g_trees_df$polyrhiza)),
  H_annuus=dist_out_g_trees_df$annuus/mean(na.omit(dist_out_g_trees_df$annuus))
)

### Write a CSV of of the distances
#Add subtree, Athal, Osat info
Dist4write_df_norm<-cbind(data.frame(Subtree=c("RNApol", "ACCD", "CLPP1", "MATK", "PHOTOSYNTH", "RIBO", "YCF", unlist(Subtree_file)), 
                 Athal_seq=unlist(Athal_names),
                 Osat_seq=unlist(Osat_names)),
                 dist_out_norm_df)

Dist4write_df<-cbind(data.frame(Subtree=c("RNApol", "ACCD", "CLPP1", "MATK", "PHOTOSYNTH", "RIBO", "YCF", unlist(Subtree_file)), 
                                     Athal_seq=unlist(Athal_names),
                                     Osat_seq=unlist(Osat_names)),
                     dist_out_g_trees_df)


## Write the csv
write.csv(Dist4write_df_norm, file = "OUT_CSVs/ERC_Distance_data_norm.csv")
write.csv(Dist4write_df, file = "OUT_CSVs/ERC_Distance_data_raw.csv")

#####################################################
##########            DO ERC            #############
#####################################################

#If nessesary, read dist_norm_df data back in #NOTE: need to remove
#dist_out_norm_df<-read.csv(file = "OUT_CSVs/ERC_Distance_data_norm.csv")

#make matrix to store results about spead and p-values
spread_mat<-matrix(, ncol = 6, nrow = length(plastid_trees))

#Loop through each of the plastid trees to perform ERC analyses using each of the plast trees
for(b in 1:length(plastid_trees)){ 

  #Make blank matrix
  CxA_out<-matrix(, ncol = 8, nrow = nrow(dist_out_norm_df))
  
  #Loop through all distance measures
  for(e in 1:nrow(dist_out_norm_df)){

    #Get a dataframe that inculdes only the test genes 
    test_df<-data.frame(Branch=names(dist_out_norm_df), GeneA=as.numeric(paste(dist_out_norm_df[b,])), GeneB=as.numeric(paste(dist_out_norm_df[e,])))
   
    #remove any rows with NA
    test_df<-subset(test_df, test_df$GeneA!="NaN" & test_df$GeneB!="NaN")
  
    #Get p-values
    #pearson
    P_val_pearson<-cor.test(x = test_df$GeneA, y = test_df$GeneB, method = "pearson")$p.value
    #spearman
    P_val_spearman<-cor.test(x = test_df$GeneA, y = test_df$GeneB, method = "spearman")$p.value
    
    #Get r-squared
    r2<-summary(lm(test_df$GeneA~test_df$GeneB))$r.squared
    #get slope of trend line
    slope<-lm(test_df$GeneA~test_df$GeneB)$coefficients[2]
    
    #Fill in the output
    CxA_out[e,]<-c(as.numeric(paste(e)), as.numeric(paste(slope)), as.numeric(paste(r2)), as.numeric(paste(P_val_pearson)), as.numeric(paste(P_val_spearman)), paste(Athal_names[e]), paste(Osat_names[e]), paste(Subtree_file[e]))
  
   #  ### Plot the ERC dotplot
   #  #Note the print() below is required when I run this is the loop
   #  pdf(file = paste("PDFs_of_ERCs/",plast_names[b], "_", c(plast_names, sprintf("%04d", final_keeper_subtrees))[e], ".pdf", sep = ""), width = 4, height = 4)
   # 
   #  print(
   #    ggplot(test_df, aes(y=test_df$GeneB, x=test_df$GeneA, color=paste(test_df$Branch)))+
   #    geom_point()+
   #    geom_abline(intercept = lm(test_df$GeneB~test_df$GeneA)$coefficients[1], slope = lm(test_df$GeneB~test_df$GeneA)$coefficients[2])+
   #    xlab(plast_names[b])+
   #    ylab(substr(paste(unlist(Athal_names[e])), 29, 37))+
   #    ggtitle(paste(
   #      paste("Gene ID: ",substr(paste(unlist(Athal_names[e])), 29, 37), "\n", sep = ""),
   #      paste("Pearson p-val: ",formatC(P_val_pearson), "\n", sep = ""),
   #      paste("Spearman p-val: ",formatC(P_val_spearman), "\n", sep = ""),
   #      paste("r2: ", formatC(r2),  sep = "")
   #    ))+
   #    labs(color="Branch")+
   #    theme(plot.title = element_text(size=8, face="italic"), legend.text = element_text(size=6, face="italic"), legend.key.size = unit(0.3, "cm"))
   # )#End print
   # 
   # dev.off()
    
    } #End CvsA loop
  
  #####################################################
  ##########       Clean ERC results      #############
  #####################################################
  
  #Convert to dataframe
  CxA_out_df<-data.frame(Number=as.numeric(paste(CxA_out[,1])), 
                         Slope=as.numeric(paste(CxA_out[,2])),
                         Slope_OG=as.numeric(paste(CxA_out[,2])),
                         R_squared=as.numeric(paste(CxA_out[,3])), 
                         P_val_pearson=as.numeric(paste(CxA_out[,4])),
                         P_val_pearson_OG=as.numeric(paste(CxA_out[,4])), 
                         P_val_spearman=as.numeric(paste(CxA_out[,5])),
                         P_val_spearman_OG=as.numeric(paste(CxA_out[,5])),
                         Athal_seq=paste(CxA_out[,6]),
                         Osat_seq=paste(CxA_out[,7]),
                         Ingroup_monophyletic=c(sample(c("NA"),size = length(plastid_trees), replace = TRUE), unlist(Monophyly_test)),
                         Subtree_file=paste(CxA_out[,8]), 
                         CLP_ID=paste("Other"), 
                         CPMTloc=paste("Not_CPMTloc"))
  
  #Remove the self-by-self
  CxA_out_df<-CxA_out_df[-b,]
  
  #Convert p=0 to p = 2.2e-16
  CxA_out_df$P_val_spearman[which(CxA_out_df$P_val_spearman == 0)]<-as.numeric(2.2e-16)
  
  #We aren't interested in negative correlations so we change the slope and p-values for negative slopes to 0 and 1, respectively
  CxA_out_df$Slope[which(CxA_out_df$Slope < 0)]<-as.numeric(0)
  CxA_out_df$P_val_spearman[which(CxA_out_df$Slope_OG < 0)]<-as.numeric(1)
  CxA_out_df$P_val_pearson[which(CxA_out_df$Slope_OG < 0)]<-as.numeric(1)
  
  #Add multi-test correction
  CxA_out_df_cor<-cbind(CxA_out_df,
                    data.frame(P_value_pearson_adj=p.adjust(p=as.numeric(paste(CxA_out_df$P_val_pearson)), method = "fdr")),
                    data.frame(P_value_spearman_adj=p.adjust(p=as.numeric(paste(CxA_out_df$P_val_spearman)), method = "fdr"))
  )
  
  ###Add functional info
  CLPRs<-c("AT1G49970", "AT1G12410", "AT1G09130", "AT4G17040")
  CLPPs<-c("ATCG00670", "AT1G66670", "AT5G45390", "AT1G02560", "AT1G11750")
  CLP_others<-c("AT5G51070", "AT5G50920", "AT3G48870", "AT4G25370", "AT4G12060", "AT1G68660", "AT2G03390")
  #See Yu et al 2008 for Pseudouridine synthase info
  PSs<-c(
    "AT1G09800","AT1G20370","AT1G34150","AT1G56345","AT1G76050","AT1G76120","AT1G78910","AT2G30320","AT3G04820","AT3G06950","AT3G19440","AT3G57150","AT4G21770","AT5G14460","AT5G35400","AT5G51140"
  )
  
  #Get Cymira csv to add to plots
  cymira<-read.csv(file = "Cymira_download_v1.csv", header = TRUE)
  
  mtRibo<-paste(subset(cymira, cymira$CyMIRA.Interaction.Category == "Mitoribosome")$AGI.Identifier)
  cpRibo<-paste(subset(cymira, cymira$CyMIRA.Interaction.Category == "Chlororibosome")$AGI.Identifier)
  PPRs<-paste(subset(cymira, cymira$CyMIRA.Interaction.Category == "PPR")$AGI.Identifier)
  cpLoc<-paste(subset(cymira, cymira$CyMIRA.targeting == "Plastid")$AGI.Identifier)
  mtLoc<-paste(subset(cymira, cymira$CyMIRA.targeting == "Mitochondria")$AGI.Identifier)
  dualLoc<-paste(subset(cymira, cymira$CyMIRA.targeting == "Dual")$AGI.Identifier)
  
  #lists for results
  CLPRs_found<-list()
  CLPPs_found<-list()
  CLP_others_found<-list()
  PSs_found<-list()
  cpRibo_found<-list()
  mtRibo_found<-list()
  ppr_found<-list()
  cpLoc_found<-list()
  mtLoc_found<-list()
  dualLoc_found<-list()
  #lists for the row numbers that are hits
  CLPRs_found_n<-list()
  CLPPs_found_n<-list()
  CLP_others_found_n<-list()
  PSs_found_n<-list()
  cpRibo_found_n<-list()
  mtRibo_found_n<-list()
  ppr_found_n<-list()
  cpLoc_found_n<-list()
  mtLoc_found_n<-list()
  dualLoc_found_n<-list()
  
  #loops
  #CLPRs
  for(r in 1:length(CLPRs)){
    hit_r<-grep(CLPRs[r], CxA_out_df$Athal_seq)
    if(length(hit_r) > 0){
      CLPRs_found<-c(CLPRs_found, paste(CxA_out_df$Athal_seq[hit_r]))
      CLPRs_found_n<-c(CLPRs_found_n, hit_r)
    }
  }
  
  #CLPPs
  for(p in 1:length(CLPPs)){
    hit_p<-grep(CLPPs[p], CxA_out_df$Athal_seq)
    if(length(hit_p) > 0){
      CLPPs_found<-c(CLPPs_found, paste(CxA_out_df$Athal_seq[hit_p]))
      CLPPs_found_n<-c(CLPPs_found_n, hit_p)
    }
  }
  
  #Other CLPs
  for(o in 1:length(CLP_others)){
    hit_o<-grep(CLP_others[o], CxA_out_df$Athal_seq)
    if(length(hit_o) > 0){
      CLP_others_found<-c(CLP_others_found, paste(CxA_out_df$Athal_seq[hit_o]))
      CLP_others_found_n<-c(CLP_others_found_n, hit_o)
    }
  }
  
  #PSs
  for(s in 1:length(PSs)){
    hit_s<-grep(PSs[s], CxA_out_df$Athal_seq)
    if(length(hit_s) > 0){
      PSs_found<-c(PSs_found, paste(CxA_out_df$Athal_seq[hit_s]))
      PSs_found_n<-c(PSs_found_n, hit_s)
    }
  }
  
  #cpRibosome
  for(cp in 1:length(cpRibo)){
    hit_cp<-grep(cpRibo[cp], CxA_out_df$Athal_seq)
    if(length(hit_cp) > 0){
      cpRibo_found<-c(cpRibo_found, paste(CxA_out_df$Athal_seq[hit_cp]))
      cpRibo_found_n<-c(cpRibo_found_n, hit_cp)
    }
  }
  
  #mtRibosome
  for(mt in 1:length(mtRibo)){
    hit_mt<-grep(mtRibo[mt], CxA_out_df$Athal_seq)
    if(length(hit_mt) > 0){
      mtRibo_found<-c(mtRibo_found, paste(CxA_out_df$Athal_seq[hit_mt]))
      mtRibo_found_n<-c(mtRibo_found_n, hit_mt)
    }
  }
  
  #PPR
  for(ppr in 1:length(PPRs)){
    hit_ppr<-grep(PPRs[ppr], CxA_out_df$Athal_seq)
    if(length(hit_ppr) > 0){
      ppr_found<-c(ppr_found, paste(CxA_out_df$Athal_seq[hit_ppr]))
      ppr_found_n<-c(ppr_found_n, hit_ppr)
    }
  }
  
  #Assign CLP status
  levels(CxA_out_df_cor$CLP_ID)<-c("Other", "CLPR", "CLPP", "CLP_other", "PS", "cpRibosome", "mtRibosome", "PPR")
  CxA_out_df_cor$CLP_ID[c(as.numeric(paste(CLPRs_found_n)))]<-"CLPR"
  CxA_out_df_cor$CLP_ID[as.numeric(paste(CLPPs_found_n))]<-"CLPP"
  CxA_out_df_cor$CLP_ID[as.numeric(paste(CLP_others_found_n))]<-"CLP_other"
  CxA_out_df_cor$CLP_ID[as.numeric(paste(PSs_found_n))]<-"PS"
  CxA_out_df_cor$CLP_ID[as.numeric(paste(cpRibo_found_n))]<-"cpRibosome"
  CxA_out_df_cor$CLP_ID[as.numeric(paste(mtRibo_found_n))]<-"mtRibosome"
  CxA_out_df_cor$CLP_ID[as.numeric(paste(ppr_found_n))]<-"PPR"
  
  #Loops for localization
  #cpLoc
  for(cpL in 1:length(cpLoc)){
    hit_cpL<-grep(cpLoc[cpL], CxA_out_df_cor$Athal_seq)
    if(length(hit_cpL) > 0){
      cpLoc_found<-c(cpLoc_found, paste(CxA_out_df_cor$Athal_seq[hit_cpL]))
      cpLoc_found_n<-c(cpLoc_found_n, hit_cpL)
    }
  }
  
  #mtLoc
  for(mtL in 1:length(mtLoc)){
    hit_mtL<-grep(mtLoc[mtL], CxA_out_df_cor$Athal_seq)
    if(length(hit_mtL) > 0){
      mtLoc_found<-c(mtLoc_found, paste(CxA_out_df_cor$Athal_seq[hit_mtL]))
      mtLoc_found_n<-c(mtLoc_found_n, hit_mtL)
    }
  }
  
  #dualLoc
  for(dualL in 1:length(dualLoc)){
    hit_dualL<-grep(dualLoc[dualL], CxA_out_df_cor$Athal_seq)
    if(length(hit_dualL) > 0){
      dualLoc_found<-c(dualLoc_found, paste(CxA_out_df_cor$Athal_seq[hit_dualL]))
      dualLoc_found_n<-c(dualLoc_found_n, hit_dualL)
    }
  }
  
  #Assign Localization status
  levels(CxA_out_df_cor$CPMTloc)<-c("Not_CPMTloc", "cpLoc", "mtLoc", "dual_Loc")
  CxA_out_df_cor$CPMTloc[as.numeric(paste(cpLoc_found_n))]<-"cpLoc"
  CxA_out_df_cor$CPMTloc[as.numeric(paste(mtLoc_found_n))]<-"mtLoc"
  CxA_out_df_cor$CPMTloc[as.numeric(paste(dualLoc_found_n))]<-"dual_Loc"
  
  #Reorder by pvalue
  CxA_out_reorder_df<-CxA_out_df_cor[order(CxA_out_df_cor$P_val_pearson), ]
  
  #Write csv
  write.csv(CxA_out_reorder_df, file = paste("OUT_CSVs/",
                                             substr(
                                               plastid_trees[[b]]$tip.label[grep("thaliana", plastid_trees[[b]]$tip.label)], 22, nchar(plastid_trees[[b]]$tip.label[grep("thaliana", plastid_trees[[b]]$tip.label)])
                                             ), "_ERC_out.csv",
                                             sep = ""), quote = FALSE)
}#End plastid genes loop

