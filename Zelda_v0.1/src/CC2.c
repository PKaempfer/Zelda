/*
 * CC2.c
 *
 *  Created on: Mar 13, 2017
 *      Author: lkaempfpp
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "CC.h"
#include "ConsensusCaller.h"
#include "kmer.h"
#include "DBGraph_scaffold.h"

#define maxAltLen 10
//#define mask_N

long double POG_variant_qualaty(int n,int k){
	long double right = 0.995;
	long double wrong = 0.005;
    long long ans=1;
    k=k>n-k?n-k:k;
    int j=1;
    for(;j<=k;j++,n--)
    {
        if(n%j==0)
        {
            ans*=n/j;
        }else
        if(ans%j==0)
        {
            ans=ans/j*n;
        }else
        {
            ans=(ans*n)/j;
        }
    }
    return (-10 * log10((ans*(pow(wrong,k)*pow(right,n-k)))));
}


static inline unsigned char POG_makeCigar(char* cigar, char* ref, char* alt){
	unsigned char numM = 0;
	unsigned char numX = 0;
	unsigned char totX = 0;
	uint16_t rl = strlen(ref);
	uint16_t al = strlen(alt);
	cigar[0]='\0';
	// SNP
	if(al==1 && rl ==1){
		sprintf(cigar,"1X");
		return 0;
	}
	else{

	}
	int min = _min(rl,al);
	int i;
	for(i=0;i<min;i++){
		if(ref[i] == alt[i]){
			if(numM) numM++;
			else{
				if(numX){
					sprintf(cigar,"%s%iX",cigar,(int)numX);
					numX=0;
				}
				numM++;
			}
		}
		else{
			if(numX) numX++;
			else{
				if(numM){
					sprintf(cigar,"%s%iM",cigar,(int)numM);
					numM=0;
				}
				numX++;
			}
			totX++;
		}
	}
	if(numM) sprintf(cigar,"%s%iM",cigar,(int)numM);
	if(numX) sprintf(cigar,"%s%iX",cigar,(int)numX);
	if(rl == al) return 1;
	else{
		if(al > rl){
			sprintf(cigar,"%s%iI",cigar,al-rl);
			if(totX) return 4;
			else return 2;
		}
		else{
			sprintf(cigar,"%s%iD",cigar,rl-al);
			if(totX) return 4;
			else return 3;
		}
	}
	return 0;
}

void POG_deleteVariant(struct POG* pog){
	int i;
	struct Variation* var;
	for(i=0;i<pog->contigNum;i++){
		var = pog->contig[i].var;
		while(var){
			pog->contig[i].var = var->next;
			free(var->altSeq);
			free(var->refSeq);
			free(var);
			var = pog->contig[i].var;
		}
	}
}

void POG_reportVariant(struct POG* pog, char* vcfFile, char* ref){
	// Write Alternative to VCF
	FILE* vcf = fopen(vcfFile,"w");

	// ToDo: Write Header
    time_t current_time;
    char* cigar = (char*)malloc(2*maxAltLen);
    char* c_time_string;
    /* Obtain current time. */
    current_time = time(NULL);

    if (current_time == ((time_t)-1))
    {
        (void) fprintf(stderr, "Failure to obtain the current time.\n");
        exit(EXIT_FAILURE);
    }

    /* Convert to local time format. */
    c_time_string = ctime(&current_time);

    if (c_time_string == NULL)
    {
        (void) fprintf(stderr, "Failure to convert the current time.\n");
        exit(EXIT_FAILURE);
    }
    c_time_string[strlen(c_time_string)-1] = '\0';
    fprintf(vcf,"##fileformat=VCFv4.1\n");
    fprintf(vcf,"##fileDate=\"%s\"\n",c_time_string);
    fprintf(vcf,"##source=%s\n",version);
    fprintf(vcf,"##reference=%s\n",ref);
    fprintf(vcf,"##phasing=none\n");
    fprintf(vcf,"##filter=\"QUAL > 20\"\n");
    fprintf(vcf,"##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
    fprintf(vcf,"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">\n");
//    fprintf(vcf,"##INFO=<ID=DPB,Number=1,Type=Float,Description=\"Total read depth per bp at the locus; bases in reads overlapping / bases in haplotype\">\n");
//    fprintf(vcf,"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">\n");
//    fprintf(vcf,"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n");
//    fprintf(vcf,"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1]\">\n");
    fprintf(vcf,"##INFO=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count, with partial observations recorded fractionally\">\n");
    fprintf(vcf,"##INFO=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observations, with partial observations recorded fractionally\">\n");
//    fprintf(vcf,"##INFO=<ID=PRO,Number=1,Type=Float,Description=\"Reference allele observation count, with partial observations recorded fractionally\">\n");
//    fprintf(vcf,"##INFO=<ID=PAO,Number=A,Type=Float,Description=\"Alternate allele observations, with partial observations recorded fractionally\">\n");
//    fprintf(vcf,"##INFO=<ID=QR,Number=1,Type=Integer,Description=\"Reference allele quality sum in phred\">\n");
//    fprintf(vcf,"##INFO=<ID=QA,Number=A,Type=Integer,Description=\"Alternate allele quality sum in phred\">\n");
//    fprintf(vcf,"##INFO=<ID=PQR,Number=1,Type=Float,Description=\"Reference allele quality sum in phred for partial observations\">\n");
//    fprintf(vcf,"##INFO=<ID=PQA,Number=A,Type=Float,Description=\"Alternate allele quality sum in phred for partial observations\">\n");
//    fprintf(vcf,"##INFO=<ID=SRF,Number=1,Type=Integer,Description=\"Number of reference observations on the forward strand\">\n");
//    fprintf(vcf,"##INFO=<ID=SRR,Number=1,Type=Integer,Description=\"Number of reference observations on the reverse strand\">\n");
//    fprintf(vcf,"##INFO=<ID=SAF,Number=A,Type=Integer,Description=\"Number of alternate observations on the forward strand\">\n");
//    fprintf(vcf,"##INFO=<ID=SAR,Number=A,Type=Integer,Description=\"Number of alternate observations on the reverse strand\">\n");
//    fprintf(vcf,"##INFO=<ID=SRP,Number=1,Type=Float,Description=\"Strand balance probability for the reference allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SRF and SRR given E(SRF/SRR) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=SAP,Number=A,Type=Float,Description=\"Strand balance probability for the alternate allele: Phred-scaled upper-bounds estimate of the probability of observing the deviation between SAF and SAR given E(SAF/SAR) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=AB,Number=A,Type=Float,Description=\"Allele balance at heterozygous sites: a number between 0 and 1 representing the ratio of reads showing the reference allele to all reads, considering only reads from individuals called as heterozygous\">\n");
//    fprintf(vcf,"##INFO=<ID=ABP,Number=A,Type=Float,Description=\"Allele balance probability at heterozygous sites: Phred-scaled upper-bounds estimate of the probability of observing the deviation between ABR and ABA given E(ABR/ABA) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=RUN,Number=A,Type=Integer,Description=\"Run length: the number of consecutive repeats of the alternate allele in the reference genome\">\n");
//    fprintf(vcf,"##INFO=<ID=RPP,Number=A,Type=Float,Description=\"Read Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=RPPR,Number=1,Type=Float,Description=\"Read Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between RPL and RPR given E(RPL/RPR) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=RPL,Number=A,Type=Float,Description=\"Reads Placed Left: number of reads supporting the alternate balanced to the left (5') of the alternate allele\">\n");
//    fprintf(vcf,"##INFO=<ID=RPR,Number=A,Type=Float,Description=\"Reads Placed Right: number of reads supporting the alternate balanced to the right (3') of the alternate allele\">\n");
//    fprintf(vcf,"##INFO=<ID=EPP,Number=A,Type=Float,Description=\"End Placement Probability: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=EPPR,Number=1,Type=Float,Description=\"End Placement Probability for reference observations: Phred-scaled upper-bounds estimate of the probability of observing the deviation between EL and ER given E(EL/ER) ~ 0.5, derived using Hoeffding's inequality\">\n");
//    fprintf(vcf,"##INFO=<ID=DPRA,Number=A,Type=Float,Description=\"Alternate allele depth ratio.  Ratio between depth in samples with each called alternate allele and those without.\">\n");
//    fprintf(vcf,"##INFO=<ID=ODDS,Number=1,Type=Float,Description=\"The log odds ratio of the best genotype combination to the second-best.\">\n");
//    fprintf(vcf,"##INFO=<ID=GTI,Number=1,Type=Integer,Description=\"Number of genotyping iterations required to reach convergence or bailout.\">\n");
    fprintf(vcf,"##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">\n");
    fprintf(vcf,"##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"The extended CIGAR representation of each alternate allele, with the exception that '=' is replaced by 'M' to ease VCF parsing.\">\n");
//    fprintf(vcf,"##INFO=<ID=NUMALT,Number=1,Type=Integer,Description=\"Number of unique non-reference alleles in called genotypes at this position.\">\n");
//    fprintf(vcf,"##INFO=<ID=MEANALT,Number=A,Type=Float,Description=\"Mean number of unique non-reference allele observations per sample with the corresponding alternate alleles.\">\n");
    fprintf(vcf,"##INFO=<ID=LEN,Number=A,Type=Integer,Description=\"allele length\">\n");
//    fprintf(vcf,"##INFO=<ID=MQM,Number=A,Type=Float,Description=\"Mean mapping quality of observed alternate alleles\">\n");
//    fprintf(vcf,"##INFO=<ID=MQMR,Number=1,Type=Float,Description=\"Mean mapping quality of observed reference alleles\">\n");
//    fprintf(vcf,"##INFO=<ID=PAIRED,Number=A,Type=Float,Description=\"Proportion of observed alternate alleles which are supported by properly paired read fragments\">\n");
//    fprintf(vcf,"##INFO=<ID=PAIREDR,Number=1,Type=Float,Description=\"Proportion of observed reference alleles which are supported by properly paired read fragments\">\n");
    fprintf(vcf,"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
//    fprintf(vcf,"##FORMAT=<ID=GQ,Number=1,Type=Float,Description=\"Genotype Quality, the Phred-scaled marginal (or unconditional) probability of the called genotype\">\n");
//    fprintf(vcf,"##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood, log10-scaled likelihoods of the data given the called genotype for each possible genotype generated from the reference and alternate alleles given the sample ploidy\">\n");
    fprintf(vcf,"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
    fprintf(vcf,"##FORMAT=<ID=RO,Number=1,Type=Integer,Description=\"Reference allele observation count\">\n");
    fprintf(vcf,"##FORMAT=<ID=QR,Number=1,Type=Integer,Description=\"Sum of quality of the reference observations\">\n");
    fprintf(vcf,"##FORMAT=<ID=AO,Number=A,Type=Integer,Description=\"Alternate allele observation count\">\n");
    fprintf(vcf,"##FORMAT=<ID=QA,Number=A,Type=Integer,Description=\"Sum of quality of the alternate observations\">\n");
    fprintf(vcf,"#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  unknown\n");

	int i;
	long double qual;
	struct Variation* var;
	for(i=0;i<pog->contigNum;i++){
		var = pog->contig[i].var;
		while(var){
			qual = POG_variant_qualaty(var->dp,var->ao);
			if(qual>20){
				var->type = POG_makeCigar(cigar,var->refSeq,var->altSeq);
				fprintf(vcf,"%s\t",pog->contig[i].name);
				fprintf(vcf,"%i\t",var->pos);
				fprintf(vcf,".\t");
				fprintf(vcf,"%s\t",var->refSeq);
				fprintf(vcf,"%s\t",var->altSeq);
				fprintf(vcf,"%.3Lf\t",qual);
				fprintf(vcf,".\t");
				fprintf(vcf,"AO=%i;CIGAR=%s;DP=%i;LEN=%i;RO=%i;TYPE=%s\t",var->ao,cigar,var->dp,var->len,var->ro,varType[(int)var->type]);
				fprintf(vcf,"GT:DP:RO:QR:AO:QA\t");
				fprintf(vcf,"%i:%i:%i:%i:%i:%i\t",1,var->dp,var->ro,1,var->ao,1);
				fprintf(vcf,"\n");
			}
			var = var->next;
		}
	}

	free(cigar);

	fclose(vcf);
}

void POG_recMainPath(struct Letter_T* currentLetter, struct Letter_T* endLetter, int startPos, int altLen, char* altSeq, int altCov, struct POGseq* contig){
	char verbose = 0;
	if(verbose) printf("Checkpoint recMainPath");
	char* refSeq = (char*)malloc(maxAltLen*2);
	int refLen = 0;
	int refCov = 100000;
	struct LetterEdge* edge;
	struct Letter_T* startLetter = currentLetter;
	while(currentLetter != endLetter){
		if(refLen == maxAltLen){
			if(verbose) printf("No Proper Ref Path Found (Length Error)\n");
			return;
		}
		refSeq[refLen++] = currentLetter->letter;
		edge = currentLetter->right;
		if(edge){
			while(edge){
				if(edge->vFlag){
					refCov = _min(refCov,edge->counter);
					break;
				}
				edge = edge->next;
			}
			currentLetter = & Letters[edge->dest];
		}
	}
	char type = -1;
	if(altLen == 2 && refLen == 2) type = 0; 	// snp
	else if(altLen == refLen) type = 1; 		// mnp
	else if(altLen > refLen) type = 2;			// ins
	else if(altLen < refLen) type = 3;			// del

	struct Variation* var = (struct Variation*)malloc(sizeof(struct Variation));
	if(type == 0){
		var->altSeq = (char*)malloc(2);
		var->altSeq[0] = altSeq[1];
		var->altSeq[1] = '\0';
		var->refSeq = (char*)malloc(2);
		var->refSeq[0] = refSeq[1];
		var->refSeq[1] = '\0';
	}
	else{
		altSeq[altLen] = '\0';
		refSeq[refLen] = '\0';
		var->altSeq = (char*)malloc(strlen(altSeq)+1);
		var->refSeq = (char*)malloc(strlen(refSeq)+1);
		strcpy(var->altSeq,altSeq);
		strcpy(var->refSeq,refSeq);
	}
	var->dp = _min(startLetter->counter,endLetter->counter);
	var->ao = altCov;
	var->ro = refCov;
	var->pos = startPos;
	if(type == 0) var->pos++;
	var->next = NULL;
	var->type = type;
	// ToDo: rest entries down under
	var->len = strlen(var->altSeq);

	if(!contig->var){
		contig->var = var;
	}
	else contig->lastvar->next = var;
	contig->lastvar = var;
	if(verbose) printf("Write Variation:\n");
	if(verbose) printf("\t%s\t%i\t%s\t%s\tDP:%i;AO:%i;RO:%i\n",contig->name,var->pos,var->refSeq,var->altSeq,var->dp,var->ao,var->ro);

	free(refSeq);
	// walk main path to Letter is End of Alt Path
	// call poa_reportVariant()
}

/**
 * DFS-Search
 */
void POG_recVariantPath(struct Letter_T* startLetter, int startPos, int len, char* seq, struct LetterEdge* edge, int cov, struct POGseq* contig){
	// Include Path counting
	struct Letter_T* current = &Letters[edge->dest];
	if(len>20) printf("--> Very Long Variation (%i) at pos: %i\n",len,startPos);
	// IF current is Consensus Path call poa_recMainPath()
	if(current->vFlag) POG_recMainPath(startLetter,current,startPos,len,seq,cov,contig);
	// ELSE next not flagged go deeper
	else{
		edge = current->right;
		if(edge){
			seq[len] = current->letter;
			while(edge){
				if(edge->counter > 2 && len < maxAltLen){
					POG_recVariantPath(startLetter,startPos,len+1,seq,edge,_min(cov,edge->counter),contig);
				}
				edge = edge->next;
			}
		}
	}
}

void POG_variantCalling(struct POGseq* contig){
	static int varNum = 0;
//	printf("CHECKPOINT: Variation Calling!\n");
	int i=0;
	struct Letter_T* current = &Letters[contig->startLetter.dest];
	struct LetterEdge* edge;
	char* altPath = (char*)malloc(2*maxAltLen);
	contig->var = NULL;
	contig->lastvar = NULL;

	while(1){
		edge = current->right;
		if(edge){
			while(edge){
				if(!edge->vFlag && edge->counter > 2){
					altPath[0] = current->letter;
					varNum++;
//					if(varNum%100==0){
//						printf("Number of Variants on this path: %i\n",varNum);
//					}
					POG_recVariantPath(current,i,1,altPath,edge,edge->counter,contig);
				}
				edge = edge->next;
			}
			edge = current->right;
			while(edge){
				if(edge->vFlag) break;
				edge = edge->next;
			}
			current = &Letters[edge->dest];
			i++;
		}
		else{
			break;
		}
	}
	free(altPath);
}

void POG_avgCov(struct POGseq* contig){
	int i;
	uint32_t totCov = 0;
	for(i=0;i<numNodes;i++){
		totCov += Letters[i].counter;
	}
	contig->avgCov = (float)totCov/contig->length;
	printf("Contig: %s \t\t AvgCov: %0.2f\n",contig->name,contig->avgCov);
}

static inline void resetLetterSt(struct LetterEdge** letters){
	letters[0] = NULL;
	letters[1] = NULL;
	letters[2] = NULL;
	letters[3] = NULL;
}

void POG_alignConsensus(struct POGseq* contig, char minverbose){
//	printf("CHECKPOINT: PO-MSA to Contig\n");
	char verboseImp = 0;
	char verbose = 0;
//	if(strcmp(contig->name,"Scaffold_19_90086_186097_len:")==0) verbose = 1;
	char* seq = (char*)malloc(contig->length + 10000);
	int i=0;

	struct Letter_T* current = &Letters[contig->startLetter.dest];
	struct LetterEdge* edge;
	struct LetterEdge* bestedge;
	struct Letter_T* ring;
	struct Letter_T* bestRing = NULL;

	if(minverbose) printf("Contig: %s%i\n",contig->name,contig->length);

	struct LetterEdge* lettersSt[4];

	uint32_t backb_pos = 0;

	while(1){
//		printf("Seqpos: %i",i);
#ifdef mask_N
		if(current->counter > contig->avgCov*1.3) seq[i++] = 'N';
		else{
#endif
			if(current->letter != 'A' && current->letter != 'C' && current->letter != 'T' && current->letter != 'G' && current->letter != 'N'){
				printf("Set NON-ACGT Letter : (%c)\n",current->letter);
			}
			if(current->counter < 5) seq[i++] = current->letter+32;
			else seq[i++] = current->letter;
#ifdef mask_N
		}
#endif
//		if(current->vFlag) printf("Flag was set\n");
		current->vFlag = 1;
		edge = current->right;
		resetLetterSt(lettersSt);
		if(edge){
			bestedge = edge;
			while(edge){
				lettersSt[codes[(int)Letters[edge->dest].letter]] = edge;
				if(bestedge->counter < edge->counter){
					bestedge = edge;
				}
				edge = edge->next;
			}
			edge = bestedge;
			if(bestedge->dest < contig->length) backb_pos = bestedge->dest;
			bestRing = &Letters[bestedge->dest];
			if(Letters[bestedge->dest].align_ring){
				ring = Letters[bestedge->dest].align_ring;
				if(verbose) printf("Pos: %i -> BestRing Lettter: %c (code: %i)-> Ring Counter: %i\n",i-1,bestRing->letter,(int)bestRing->letter,bestRing->counter);
				while(ring != &Letters[bestedge->dest]){
					if(verbose) printf("\tPos: %i -> Ring Lettter: %c (code: %i)-> Ring Counter: %i\n",i-1,ring->letter,codes[(int)ring->letter],ring->counter);
					if(lettersSt[codes[(int)ring->letter]] && ring->counter > bestRing->counter && ring == &Letters[lettersSt[codes[(int)ring->letter]]->dest]){
						bestRing = ring;
						edge = lettersSt[codes[(int)ring->letter]];
					}
					ring = ring->align_ring;
				}
				if(verbose) printf("Pos: %i -> BestRing Lettter: %c (code: %i)-> Ring Counter: %i\n",i-1,bestRing->letter,(int)bestRing->letter,bestRing->counter);
			}
			edge->vFlag = 1;
			current = bestRing;
		}
		else{
			current->vFlag = 1;
			if(verbose) printf("Break at POS: %i\n",i);
			if(verbose) printf("Search for alignment ring\n");
//			bestRing = current;
//			current = current->align_ring;
//			while(current && current != bestRing){
//				if(current->right){
//					edge = current->right;
//					break;
//				}
//				current = current->align_ring;
//			}
			if(!edge){
				backb_pos++;
				if(backb_pos < (contig->length/4)*3){
					// Todo: Important to deal better with this case: Insert N's could be a proper solution.
					if(verboseImp) printf("No edge found, break CC\n");
					if(verboseImp) printf("Jump to back to Backbone (%i)\n",backb_pos);
					current = &Letters[backb_pos];
					bestRing = NULL;
					continue;
				}
				break;
			}
			else printf("New edge found, continue CC\n");

		}
		bestRing = NULL;
	}

	seq[i]='\0';
//	printf(">Correct_%s\n",contig->name);
//	int k;
//	for(k=0;k<i;k+=80){
//		printf("%.80s\n",&seq[k]);
//	}
//	printf("\n");
	contig->length = strlen(seq);
	contig->sequence = (char*)malloc(sizeof(char)*(strlen(seq)+1));
	strcpy(contig->sequence,seq);
	free(seq);
//	contig->sequence = seq;
	sprintf(contig->name,"%s%i",contig->name,contig->length);
//	struct rusage r_usage;
//	getrusage(RUSAGE_SELF,&r_usage);
//	printf("Memory usage: %ld bytes\n",r_usage.ru_maxrss);
	POG_variantCalling(contig);
//	POG_avgCov(contig);
	sprintf(contig->name,"%s_avgCov:%.2f",contig->name,contig->avgCov);
//	if(verbose) exit(1);
}


/**
 * Function prints the alignment matrix according to the poa graph. Following lines could form alignment rings and may not be successors.
 * @param row		Number of the letters in the read
 * @param column	Number of the letters in the part of the reference
 * @param seq		Sequence of the read in correct orientation
 */
void POG_showMatrix(int row, int column, char* seq){
	printf("CHECKPOINT: Print Matrix (row: %i, col: %i)\n",row,column);
	int i,j;
	printf("\t");
	for(i=0;i<strlen(seq);i++){
		printf("\t%c",seq[i]);
	}
	printf("\n");

	for(i=0;i<=row;i++){
		printf("\t%i",alMatrix[0][i]);
	}
	printf("\n");
	for(j=1;j<column;j++){
		printf("%c (%i)",alMatrix_Letter[j]->letter,j);
		for(i=0;i<=row;i++){
			printf("\t%i",alMatrix[j][i]);
		}
		printf("\n");
	}
}
