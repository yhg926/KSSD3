# --- 1) Compute alpha0_bar from FASTQ (handles variable lengths) ---
FASTQ=GCA_025985485.fq   # or .fq.gz (use zcat)
K=28
A0_100=16.364969175
A0_150=15.7789105059

# length histogram + weighted alpha0(L)
awk -v k="$K" -v a100="$A0_100" -v a150="$A0_150" '
  NR%4==2 { L=length($0); if (L>=k) {
      w = L - k + 1
      a0 = a100 + (L-100.0)*(a150-a100)/50.0   # linear interp 100->150
      W  += w
      A  += w * a0
  } }
  END { printf("%.10f\n", A/W) }
' "$FASTQ" > .alpha0_bar.txt

ALPHA0_BAR=$(cat .alpha0_bar.txt)
BETA=6.316020058

# --- 2) Get H1, H2, I from the KSSD sketch histogram ---
SKETCH=GCA_025985485_testsk    # your sketch folder
kssd3 set --psketch "$SKETCH" | cut -f2- | datamash -s -g 2 countunique 1 \
| awk -v a0="$ALPHA0_BAR" -v b="$BETA" '
    { H[$1]=$2; I += $1*$2 }
    END{
      H1 = (H[1]+0); H2 = (H[2]+0);
      s  = (H1 - 2*H2)/I;
      r  = s / (a0 - b*s);
      printf "alpha0_bar=%.6f  H1=%d  H2=%d  I=%d  r=%.6f (%.2f%%)\n",
             a0, H1,H2,I, r, 100*r;
    }'

