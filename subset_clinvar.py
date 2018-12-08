#!/usr/bin/python
import sys
import gzip

outfile=open(sys.argv[1], "w")

clinvar=gzip.open("/clinvar/output/b37/single/clinvar_alleles_with_gnomad_exomes.single.b37.tsv.gz", 'rb')
header=clinvar.readline().rstrip('\n').split('\t')
keep_cols=["chrom","pos","ref","alt","allele_id","symbol","hgvs_c","hgvs_p","molecular_consequence","clinical_significance","review_status","all_traits","AC", "AN", "AF","Hom","AC_AFR","AC_AMR","AC_ASJ","AC_EAS","AC_SAS","AC_FIN","AC_NFE","AN_AFR","AN_AMR","AN_ASJ","AN_EAS","AN_SAS","AN_FIN","AN_NFE","Hom_AFR","Hom_AMR","Hom_ASJ","Hom_EAS","Hom_SAS","Hom_FIN","Hom_NFE","AF_POPMAX"]
clnsig_col=header.index("clinical_significance")
cons_col=header.index("molecular_consequence")
indexes = [index for index in range(len(header)) if header[index] in keep_cols]
clinvar.close()


clinvar=gzip.open("/clinvar/output/b37/single/clinvar_alleles_with_gnomad_exomes.single.b37.tsv.gz", 'rb')
for line_c1 in clinvar:
        line_c=line_c1.rstrip('\n').split('\t')
        if line_c[0]!="chrom":
                clnsig=line_c[clnsig_col]
                cons=line_c[cons_col]
                if ("athog" in clnsig) and ("risk" not in clnsig) and ("drug" not in clnsig) and ("enign" not in clnsig) and ("flict" not in clnsig) and ("certain" not in clnsig):
                        if "hello" not in cons:
                                if str(line_c[0]) not in ["X", "MT", "Y"]:
                                        outfile.write('\t'.join([line_c[i] for i in indexes])+"\n")
        else:
                outfile.write('\t'.join([line_c[i] for i in indexes])+"\n")
clinvar.close()
outfile.close()
