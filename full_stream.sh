#!/bin/bash

# updated 2023/04/05 by JCT

ID=$1
corepath=$"/gpfs3/well/win-biobank/projects/imaging/data/data3/subjectsAll/""$ID""/"
homepath=$"/well/webb/users/vox025/"

echo $ID start full stream `date`

mkdir -p "$homepath""$ID""_analysis/masks/"
mkdir -p "$homepath""$ID""_analysis/preprocessing/"
mkdir -p "$homepath""$ID""_analysis/alff/alff0/"
mkdir -p "$homepath""$ID""_analysis/alff/alff1/"
mkdir -p "$homepath""$ID""_analysis/alff/alff2/"
mkdir -p "$homepath""$ID""_analysis/alff/alff3/"
mkdir -p "$homepath""$ID""_analysis/alff/alff4/"

###################################################################################################################################################################################################

echo ----- $ID start file check `date`

cd "$homepath""$ID""_analysis/preprocessing/"

file_check="$corepath""/T1/T1.nii.gz"
if [ -f "$file_check" ]; then
	T1=1
else
	T1=0
fi

file_check="$corepath""/T2_FLAIR/T2_FLAIR.nii.gz"
if [ -f "$file_check" ]; then
	T2_FLAIR=1
else
	T2_FLAIR=0
fi

file_check="$corepath""/fMRI/rfMRI.nii.gz"
if [ -f "$file_check" ]; then
	rfMRI=1
else
	rfMRI=0
fi

file_check="$corepath""/SWI/SWI_TOTAL_MAG_to_T1.nii.gz"
if [ -f "$file_check" ]; then
	SWI=1
else
	SWI=0
fi

cat >> ./"$ID"_file_check.txt << EOL
$1 $T1 $T2_FLAIR $rfMRI $SWI
EOL

cat ./"$ID"_file_check.txt >> "$homepath"final_outputs/all_file_check.txt

echo -------- $ID complete file check `date`

###################################################################################################################################################################################################

echo ----- $ID start vascular masking `date`

cd "$homepath""$ID""_analysis/masks/"

# Create warp for mask containing sag sinus and pulsatile mask from standard space to T1 high res space (inverting warp in other direction in Biobank folder)
invwarp --ref="$corepath""/T1/T1.nii.gz"  --warp="$corepath""/T1/transforms/T1_to_MNI_warp_coef.nii.gz" --out="MNI_to_T1_warp.nii.gz"
applywarp --ref="$corepath""/T1/T1.nii.gz" --in="$homepath""/MNI_masks/MNI_sagsinus.nii" --warp="MNI_to_T1_warp.nii.gz" --out="T1_sagsinus.nii.gz" --interp=nn
applywarp --ref="$corepath""/T1/T1.nii.gz" --in="$homepath""/MNI_masks/pulsatility_mask_comb.nii.gz" --warp="MNI_to_T1_warp.nii.gz" --out="pulsatility_mask_comb.nii.gz" --interp=nn

# Split combined mask in structural space into a sagsinus mask and a 'Circle of Willis' mask
fslmaths pulsatility_mask_comb.nii.gz -thr 1.6 -bin pulsatility_mask.nii.gz
fslmaths pulsatility_mask_comb.nii.gz -thr 0.3 -uthr 1 -bin T1_sagsinus.nii.gz

#Expand sagsinus mask to ensure good coverage in individual
fslmaths T1_sagsinus.nii.gz -dilD T1_sagsinus_mod.nii.gz

#Refine sag sinus mask by high signal on SWI
fslmaths "$corepath""/SWI/SWI_TOTAL_MAG_to_T1.nii.gz" -mas T1_sagsinus_mod.nii.gz SWI_sag_masked.nii.gz
thresh=$(fslstats SWI_sag_masked.nii.gz -P 98)
thresh=$(bc <<< "$thresh*0.6")
fslmaths SWI_sag_masked.nii.gz -thr $thresh -bin SWI_sag_masked.nii.gz

#Refine sagsinus mask by T2
thresh2=$(fslstats "$corepath""/T2_FLAIR/T2_FLAIR.nii.gz" -P 98)
thresh2=$(bc <<< "$thresh2*0.1")
fslmaths "$corepath""/T2_FLAIR/T2_FLAIR.nii.gz" -thr $thresh2 -binv T2_mask.nii.gz
fslmaths SWI_sag_masked.nii.gz -mas T2_mask.nii.gz SWI_sag_masked.nii.gz
fslmaths SWI_sag_masked.nii.gz -dilD SWI_sag_masked.nii.gz
fslmaths SWI_sag_masked.nii.gz -ero T1_sagsinus_mod.nii.gz

#Pick out top few clusters to get rid of the small, non-connected cluster of voxels
cluster -i T1_sagsinus_mod.nii.gz -t 0.5 -o cluster_index.nii.gz
clustmax=$(fslstats cluster_index.nii.gz -P 100)
fslmaths -dt int cluster_index.nii.gz -thr $clustmax -uthr $clustmax -bin T1_sagsinus_mod.nii.gz

#Derive vascular mask by same method
fslmaths pulsatility_mask.nii.gz -nan pulsatility_mask.nii.gz
fslmaths pulsatility_mask.nii.gz -dilD pulsatility_mask_mod.nii.gz
fslmaths "$corepath""/SWI/SWI_TOTAL_MAG_to_T1.nii.gz" -mas pulsatility_mask_mod.nii.gz SWI_pulse_masked.nii.gz

thresh2=$(fslstats "$corepath""/T2_FLAIR/T2_FLAIR.nii.gz" -P 95)
thresh2=$(bc <<< "$thresh2*0.4")
fslmaths "$corepath""/T2_FLAIR/T2_FLAIR.nii.gz" -thr $thresh2 -binv T2_mask.nii.gz
fslmaths SWI_pulse_masked.nii.gz -mas T2_mask.nii.gz SWI_pulse_masked_T2.nii.gz

thresh=$(fslstats SWI_pulse_masked.nii.gz -P 75)
thresh=$(bc <<< "$thresh*0.8")

#Split extracted pulsatility mask into vessels and CSF around vessels
fslmaths SWI_pulse_masked_T2.nii.gz -thr $thresh -bin SWI_pulse_masked_tight.nii.gz
fslmaths SWI_pulse_masked_tight.nii.gz -dilD pulsatility_mask_mod.nii.gz
fslmaths SWI_pulse_masked_T2.nii.gz -mas SWI_pulse_masked_tight.nii.gz tight_value.nii.gz
fslmaths SWI_pulse_masked_T2.nii.gz -sub tight_value.nii.gz SWI_pulse_masked_T2_CSF.nii.gz
fslmaths SWI_pulse_masked_T2.nii.gz -bin SWI_pulse_masked_T2.nii.gz
fslmaths SWI_pulse_masked_T2_CSF.nii.gz -bin SWI_pulse_masked_T2_CSF.nii.gz

# Create masks in functional space
convert_xfm -omat T1_to_BOLD.mat -inverse "$corepath""/fMRI/rfMRI.ica/reg/example_func2highres.mat"

flirt -in T1_sagsinus_mod.nii.gz -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_sagsinus_mod.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  BOLD_sagsinus_mod.nii.gz -thr 0.3 -bin BOLD_sagsinus_mod.nii.gz

flirt -in "$corepath""/T1/T1_fast/T1_brain_seg.nii.gz"  -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_brain_seg.nii.gz -init T1_to_BOLD.mat -applyxfm
#fslmaths  ./BOLD_brain_seg -thr 0.3 -bin ./BOLD_brain_seg

# Pick out mask of deep grey matter from segmentation
fslmaths "$corepath""/T1/T1_first/T1_first_all_fast_firstseg.nii.gz" -uthr 13 "T1_deepGM.nii.gz"
fslmaths "$corepath""/T1/T1_first/T1_first_all_fast_firstseg.nii.gz" -uthr 52 "test.nii.gz"
fslmaths "test.nii.gz" -thr 49 "test.nii.gz"
fslmaths "T1_deepGM.nii.gz" -add "test.nii.gz" "T1_deepGM.nii.gz"
rm "test.nii.gz"

# flirt deep grey to functional space
flirt -in "T1_deepGM.nii.gz"  -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_first_seg.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  ./BOLD_first_seg.nii.gz -thr 0.3 -bin ./BOLD_first_seg.nii.gz

flirt -in pulsatility_mask_mod.nii.gz -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_pulsatility_mask_mod.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  BOLD_pulsatility_mask_mod.nii.gz -thr 0.5 -bin BOLD_pulsatility_mask_mod.nii.gz

cluster -i BOLD_pulsatility_mask_mod -t 0.5 -o cluster_index.nii.gz
clustmax=$(fslstats cluster_index.nii.gz -P 100)
clustmin=$(bc <<< "$clustmax - 4")
fslmaths -dt int cluster_index.nii.gz -thr $clustmin -uthr $clustmax -bin BOLD_pulsatility_mask_mod.nii.gz

#Sections to copy other versions of pulsatile masks into fucntional space

#Non-modified version (ie not isolated out vessels)
flirt -in pulsatility_mask.nii.gz -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_pulsatility_mask.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  BOLD_pulsatility_mask.nii.gz -thr 0.5 -bin BOLD_pulsatility_mask.nii.gz

# Version of vessels without dilatation steps - much tighter but likely to less perfect aligned and continuous once in functional space
flirt -in SWI_pulse_masked_tight.nii.gz -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_pulse_masked_tight.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  BOLD_pulse_masked_tight.nii.gz -thr 0.5 -bin BOLD_pulse_masked_tight.nii.gz

# T2 dependent mask of pulsatile segments - is pulling out the CSF
flirt -in SWI_pulse_masked_T2.nii.gz -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_pulse_masked_T2.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  BOLD_pulse_masked_T2.nii.gz -thr 0.5 -bin BOLD_pulse_masked_T2.nii.gz

# Same as T2 but removes vessels to give the CSF parts around it
flirt -in SWI_pulse_masked_T2_CSF.nii.gz -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_pulse_masked_T2_CSF.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  BOLD_pulse_masked_T2_CSF.nii.gz -thr 0.8 -bin BOLD_pulse_masked_T2_CSF.nii.gz

echo -------- $ID complete vascular masking `date`

###################################################################################################################################################################################################

echo ----- $ID start FAST masking `date`

cd "$homepath""$ID""_analysis/masks/"

flirt -in "$corepath""/T1/T1_fast/T1_brain_pve_0.nii.gz"  -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out ./BOLD_brain_CSF.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  ./BOLD_brain_CSF -thr 0.5 -bin ./BOLD_brain_CSF3

flirt -in "$corepath""/T1/T1_fast/T1_brain_pve_1.nii.gz"  -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out ./BOLD_brain_GM.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  ./BOLD_brain_GM -thr 0.6 -bin ./BOLD_brain_GM4

flirt -in "$corepath""/T1/T1_fast/T1_brain_pve_2.nii.gz"  -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out ./BOLD_brain_WM.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths  ./BOLD_brain_WM -thr 0.7 -bin ./BOLD_brain_WM5

fslmaths ./BOLD_brain_CSF3 -sub ./BOLD_brain_WM5 -sub ./BOLD_brain_GM4 -bin ./BOLD_brain_CSF_exclusive
fslmaths ./BOLD_brain_GM4 -sub ./BOLD_brain_WM5 -sub ./BOLD_brain_CSF3 -bin ./BOLD_brain_GM_exclusive
fslmaths ./BOLD_brain_WM5 -sub ./BOLD_brain_CSF3 -sub ./BOLD_brain_GM4 -bin ./BOLD_brain_WM_exclusive
fslmaths  ./BOLD_first_seg -thr 0.5 -bin ./BOLD_first_seg5

echo -------- $ID complete FAST masking `date`

###################################################################################################################################################################################################

echo ----- $ID start BIANCA masking `date`

cd "$homepath""$ID""_analysis/masks/"

flirt -in "$corepath""/T2_FLAIR/lesions/pwmh_dwmh_output/deepwm_map_10mm.nii.gz" -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_bianca_dwmh.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths BOLD_bianca_dwmh.nii.gz -thr 0.5 -bin BOLD_bianca_dwmh.nii.gz
flirt -in "$corepath""/T2_FLAIR/lesions/pwmh_dwmh_output/perivent_map_10mm.nii.gz" -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out BOLD_bianca_pwmh.nii.gz -init T1_to_BOLD.mat -applyxfm
fslmaths BOLD_bianca_pwmh.nii.gz -thr 0.5 -bin BOLD_bianca_pwmh.nii.gz
fslmaths BOLD_bianca_dwmh.nii.gz -add BOLD_bianca_pwmh.nii.gz -bin BOLD_bianca_allwmh.nii.gz
fslmaths BOLD_brain_WM_exclusive.nii.gz -sub BOLD_bianca_allwmh.nii.gz -bin BOLD_nawm.nii.gz

echo -------- $ID complete BIANCA masking `date`

###################################################################################################################################################################################################

echo ----- $ID start fMRI BET and MCFLIRT `date`

cd "$homepath""$ID""_analysis/preprocessing/"

###Motion correction
mcflirt -in "$corepath""/fMRI/rfMRI.nii.gz" -r "$corepath""/fMRI/rfMRI_SBREF.nii.gz" -o "rfMRI_mcf.nii.gz"

###Brain extraction
#bet rfMRI_mcf.nii.gz rfMRI_bet.nii.gz -F

flirt -in "$corepath""/T1/T1_brain_mask.nii.gz"  -ref "$corepath""/fMRI/rfMRI_SBREF.nii.gz"  -out ./rfMRI_bet_mask.nii.gz -init "$homepath""$ID""_analysis/masks/T1_to_BOLD.mat" -applyxfm
fslmaths  ./rfMRI_bet_mask.nii.gz -thr 0.5 -bin ./rfMRI_bet_mask.nii.gz

###Add arterial and sagittal sinus masks to bet to ensure vessels are not lost
fslmaths "rfMRI_bet_mask.nii.gz" -add "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -add "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -bin "$homepath""$ID""_analysis/masks/BOLD_bet_vasc.nii.gz"
fslmaths "rfMRI_mcf.nii.gz" -mul "$homepath""$ID""_analysis/masks/BOLD_bet_vasc.nii.gz" "rfMRI_bet_better.nii.gz"

###Create an arterial only mask (i.e. arterial subract brain)
fslmaths "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -sub "rfMRI_bet_mask.nii.gz" -bin "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" 

###Mask out deep grey from white matter
fslmaths "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive.nii.gz" -sub "$homepath""$ID""_analysis/masks/BOLD_first_seg5.nii.gz" -bin "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" 

echo -------- $ID complete fMRI BET and MCFLIRT `date`

###################################################################################################################################################################################################

echo ----- $ID start fMRI bptf `date`

cd "$homepath""$ID""_analysis/preprocessing/"

###AFNI bptf
3dTproject -stopband 0 0.009 -prefix rfMRI_bptf -input "rfMRI_bet_better.nii.gz"

##Organise ALFF files
3dAFNItoNIFTI rfMRI_bptf*
rm *.BRIK *.HEAD
gzip rfMRI_bptf*

echo -------- $ID complete fMRI bptf `date`

###################################################################################################################################################################################################

echo ----- $ID start fMRI normalisation `date`

cd "$homepath""$ID""_analysis/preprocessing/"

###Calculate time series mean, subtract mean, divide by robust range
fslmaths "rfMRI_bptf.nii.gz" -Tmean "mean_voxels.nii.gz"
fslmaths "rfMRI_bptf.nii.gz" -Tperc 5 "5_p.nii.gz"
fslmaths "rfMRI_bptf.nii.gz" -Tperc 95 "95_p.nii.gz"

fslmaths "95_p.nii.gz" -sub "5_p.nii.gz" "robust_range.nii.gz"
fslmaths "rfMRI_bptf.nii.gz" -sub "mean_voxels.nii.gz" -div "robust_range.nii.gz" "rfMRI_normalised.nii.gz"

rm 95_p* 5_p* mean_voxels* robust_range*

echo -------- $ID complete fMRI normalisation `date`

###################################################################################################################################################################################################

echo ----- $ID start masked time series `date`

cd "$homepath""$ID""_analysis/preprocessing/"

masked_ts_GM=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz") | sed -e 's/^[[:space:]]*//')
masked_ts_WM=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz") | sed -e 's/^[[:space:]]*//')
masked_ts_CSF=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz") | sed -e 's/^[[:space:]]*//')
masked_ts_arterial=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz") | sed -e 's/^[[:space:]]*//')
masked_ts_sagsinus=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz") | sed -e 's/^[[:space:]]*//')
masked_ts_allwmh=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz") | sed -e 's/^[[:space:]]*//')
masked_ts_dwmh=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz") | sed -e 's/^[[:space:]]*//')
masked_ts_pwmh=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz") | sed -e 's/^[[:space:]]*//')
masked_ts_nawm=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz") | sed -e 's/^[[:space:]]*//')
masked_ts_arterial_only=$((fslmeants -i "rfMRI_normalised.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz") | sed -e 's/^[[:space:]]*//')

cat >> ./masked_ts.txt << EOL
$1 $masked_ts_arterial $masked_ts_sagsinus $masked_ts_CSF $masked_ts_GM $masked_ts_WM $masked_ts_allwmh $masked_ts_dwmh $masked_ts_pwmh $masked_ts_nawm $masked_ts_arterial_only
EOL

cat ./masked_ts.txt >> "$homepath"final_outputs/all_masked_ts.txt

echo -------- $ID complete masked time series `date`

###################################################################################################################################################################################################

echo ----- $ID start fsl pspec `date`

cd "$homepath""$ID""_analysis/preprocessing/"

fslpspec "rfMRI_normalised.nii.gz" "rfMRI_pspec.nii.gz"

global_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose) | sed -e 's/^[[:space:]]*//')
GM_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz") | sed -e 's/^[[:space:]]*//')
WM_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz") | sed -e 's/^[[:space:]]*//')
CSF_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz") | sed -e 's/^[[:space:]]*//')
arterial_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz") | sed -e 's/^[[:space:]]*//')
sagsinus_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz") | sed -e 's/^[[:space:]]*//')
allwmh_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz") | sed -e 's/^[[:space:]]*//')
dwmh_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz") | sed -e 's/^[[:space:]]*//')
pwmh_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz") | sed -e 's/^[[:space:]]*//')
nawm_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz") | sed -e 's/^[[:space:]]*//')
arterial_only_pspec=$((fslmeants -i "rfMRI_pspec.nii.gz" --transpose -m "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz") | sed -e 's/^[[:space:]]*//')

cat >> ./"$ID"_pspec.txt << EOL
$1 $global_pspec $arterial_pspec $sagsinus_pspec $CSF_pspec $GM_pspec $WM_pspec $allwmh_pspec $dwmh_pspec $pwmh_pspec $nawm_pspec $arterial_only_pspec
EOL

cat ./"$ID"_pspec.txt >> "$homepath"final_outputs/all_pspec.txt

echo -------- $ID complete fsl pspec `date`

###################################################################################################################################################################################################

echo ----- $ID start ALFF1 `date`

cd "$homepath""$ID""_analysis/alff/alff1/"

##Create ALFF files
3dRSFC -prefix ALFF 0.05 0.1 "$homepath""$ID""_analysis/preprocessing/rfMRI_normalised.nii.gz"

##Organise ALFF files
3dAFNItoNIFTI ALFF*ALFF*; 3dAFNItoNIFTI ALFF*fALFF*; 3dAFNItoNIFTI ALFF*RSFA*; 3dAFNItoNIFTI ALFF*fRSFA*; 3dAFNItoNIFTI ALFF*mALFF*; 3dAFNItoNIFTI ALFF*mRSFA*
rm *.BRIK *.HEAD
gzip ALFF*

#Calculate ALFF stats
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -P 97.5)
ALFF1=$(fslstats ALFF_ALFF.nii.gz -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -P 97.5)
fALFF1=$(fslstats ALFF_fALFF.nii.gz -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -P 97.5)
mALFF1=$(fslstats ALFF_mALFF.nii.gz -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -P 97.5)
fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -P 97.5)
mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for sagsinus
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for CSF
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for GM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for WM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for allwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for dwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for pwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for nawm
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial only mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_ALFF1=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fALFF1=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mALFF1=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fRSFA1=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mRSFA1=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

cat >> ./"$ID"_alff1_data.txt << EOL
$1 $ALFF1 $fALFF1 $mALFF1 $fRSFA1 $mRSFA1 $arterial_ALFF1 $arterial_fALFF1 $arterial_mALFF1 $arterial_fRSFA1 $arterial_mRSFA1 $sagsinus_ALFF1 $sagsinus_fALFF1 $sagsinus_mALFF1 $sagsinus_fRSFA1 $sagsinus_mRSFA1 $CSF_ALFF1 $CSF_fALFF1 $CSF_mALFF1 $CSF_fRSFA1 $CSF_mRSFA1 $GM_ALFF1 $GM_fALFF1 $GM_mALFF1 $GM_fRSFA1 $GM_mRSFA1 $WM_ALFF1 $WM_fALFF1 $WM_mALFF1 $WM_fRSFA1 $WM_mRSFA1 $allwmh_ALFF1 $allwmh_fALFF1 $allwmh_mALFF1 $allwmh_fRSFA1 $allwmh_mRSFA1 $dwmh_ALFF1 $dwmh_fALFF1 $dwmh_mALFF1 $dwmh_fRSFA1 $dwmh_mRSFA1 $pwmh_ALFF1 $pwmh_fALFF1 $pwmh_mALFF1 $pwmh_fRSFA1 $pwmh_mRSFA1 $nawm_ALFF1 $nawm_fALFF1 $nawm_mALFF1 $nawm_fRSFA1 $nawm_mRSFA1 $arterial_only_ALFF1 $arterial_only_fALFF1 $arterial_only_mALFF1 $arterial_only_fRSFA1 $arterial_only_mRSFA1
EOL

cat ./"$ID"_alff1_data.txt >> "$homepath"final_outputs/all_alff1_data.txt

echo -------- $ID complete ALFF1 `date`

###################################################################################################################################################################################################

echo ----- $ID start ALFF2 `date`

cd "$homepath""$ID""_analysis/alff/alff2/"

##Create ALFF files
3dRSFC -prefix ALFF 0.1 0.15 "$homepath""$ID""_analysis/preprocessing/rfMRI_normalised.nii.gz"

##Organise ALFF files
3dAFNItoNIFTI ALFF*ALFF*; 3dAFNItoNIFTI ALFF*fALFF*; 3dAFNItoNIFTI ALFF*RSFA*; 3dAFNItoNIFTI ALFF*fRSFA*; 3dAFNItoNIFTI ALFF*mALFF*; 3dAFNItoNIFTI ALFF*mRSFA*
rm *.BRIK *.HEAD
gzip ALFF*

#Calculate ALFF stats
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -P 97.5)
ALFF2=$(fslstats ALFF_ALFF.nii.gz -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -P 97.5)
fALFF2=$(fslstats ALFF_fALFF.nii.gz -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -P 97.5)
mALFF2=$(fslstats ALFF_mALFF.nii.gz -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -P 97.5)
fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -P 97.5)
mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for sagsinus
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for CSF
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for GM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for WM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for allwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for dwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for pwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for nawm
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial only mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_ALFF2=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fALFF2=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mALFF2=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fRSFA2=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mRSFA2=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

cat >> ./"$ID"_alff2_data.txt << EOL
$1 $ALFF2 $fALFF2 $mALFF2 $fRSFA2 $mRSFA2 $arterial_ALFF2 $arterial_fALFF2 $arterial_mALFF2 $arterial_fRSFA2 $arterial_mRSFA2 $sagsinus_ALFF2 $sagsinus_fALFF2 $sagsinus_mALFF2 $sagsinus_fRSFA2 $sagsinus_mRSFA2 $CSF_ALFF2 $CSF_fALFF2 $CSF_mALFF2 $CSF_fRSFA2 $CSF_mRSFA2 $GM_ALFF2 $GM_fALFF2 $GM_mALFF2 $GM_fRSFA2 $GM_mRSFA2 $WM_ALFF2 $WM_fALFF2 $WM_mALFF2 $WM_fRSFA2 $WM_mRSFA2 $allwmh_ALFF2 $allwmh_fALFF2 $allwmh_mALFF2 $allwmh_fRSFA2 $allwmh_mRSFA2 $dwmh_ALFF2 $dwmh_fALFF2 $dwmh_mALFF2 $dwmh_fRSFA2 $dwmh_mRSFA2 $pwmh_ALFF2 $pwmh_fALFF2 $pwmh_mALFF2 $pwmh_fRSFA2 $pwmh_mRSFA2 $nawm_ALFF2 $nawm_fALFF2 $nawm_mALFF2 $nawm_fRSFA2 $nawm_mRSFA2 $arterial_only_ALFF2 $arterial_only_fALFF2 $arterial_only_mALFF2 $arterial_only_fRSFA2 $arterial_only_mRSFA2
EOL

cat ./"$ID"_alff2_data.txt >> "$homepath"final_outputs/all_alff2_data.txt

echo -------- $ID complete ALFF2 `date`

###################################################################################################################################################################################################

echo ----- $ID start ALFF3 `date`

cd "$homepath""$ID""_analysis/alff/alff3/"

##Create ALFF files
3dRSFC -prefix ALFF 0.05 0.15 "$homepath""$ID""_analysis/preprocessing/rfMRI_normalised.nii.gz"

##Organise ALFF files
3dAFNItoNIFTI ALFF*ALFF*; 3dAFNItoNIFTI ALFF*fALFF*; 3dAFNItoNIFTI ALFF*RSFA*; 3dAFNItoNIFTI ALFF*fRSFA*; 3dAFNItoNIFTI ALFF*mALFF*; 3dAFNItoNIFTI ALFF*mRSFA*
rm *.BRIK *.HEAD
gzip ALFF*

#Calculate ALFF stats
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -P 97.5)
ALFF3=$(fslstats ALFF_ALFF.nii.gz -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -P 97.5)
fALFF3=$(fslstats ALFF_fALFF.nii.gz -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -P 97.5)
mALFF3=$(fslstats ALFF_mALFF.nii.gz -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -P 97.5)
fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -P 97.5)
mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for sagsinus
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for CSF
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for GM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for WM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for allwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for dwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for pwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for nawm
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial only mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_ALFF3=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fALFF3=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mALFF3=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fRSFA3=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mRSFA3=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

cat >> ./"$ID"_alff3_data.txt << EOL
$1 $ALFF3 $fALFF3 $mALFF3 $fRSFA3 $mRSFA3 $arterial_ALFF3 $arterial_fALFF3 $arterial_mALFF3 $arterial_fRSFA3 $arterial_mRSFA3 $sagsinus_ALFF3 $sagsinus_fALFF3 $sagsinus_mALFF3 $sagsinus_fRSFA3 $sagsinus_mRSFA3 $CSF_ALFF3 $CSF_fALFF3 $CSF_mALFF3 $CSF_fRSFA3 $CSF_mRSFA3 $GM_ALFF3 $GM_fALFF3 $GM_mALFF3 $GM_fRSFA3 $GM_mRSFA3 $WM_ALFF3 $WM_fALFF3 $WM_mALFF3 $WM_fRSFA3 $WM_mRSFA3 $allwmh_ALFF3 $allwmh_fALFF3 $allwmh_mALFF3 $allwmh_fRSFA3 $allwmh_mRSFA3 $dwmh_ALFF3 $dwmh_fALFF3 $dwmh_mALFF3 $dwmh_fRSFA3 $dwmh_mRSFA3 $pwmh_ALFF3 $pwmh_fALFF3 $pwmh_mALFF3 $pwmh_fRSFA3 $pwmh_mRSFA3 $nawm_ALFF3 $nawm_fALFF3 $nawm_mALFF3 $nawm_fRSFA3 $nawm_mRSFA3 $arterial_only_ALFF3 $arterial_only_fALFF3 $arterial_only_mALFF3 $arterial_only_fRSFA3 $arterial_only_mRSFA3
EOL

cat ./"$ID"_alff3_data.txt >> "$homepath"final_outputs/all_alff3_data.txt

echo -------- $ID complete ALFF3 `date`

###################################################################################################################################################################################################

echo ----- $ID start ALFF4 `date`

cd "$homepath""$ID""_analysis/alff/alff4/"

##Create ALFF files
3dRSFC -prefix ALFF 0.01 0.15 "$homepath""$ID""_analysis/preprocessing/rfMRI_normalised.nii.gz"

##Organise ALFF files
3dAFNItoNIFTI ALFF*ALFF*; 3dAFNItoNIFTI ALFF*fALFF*; 3dAFNItoNIFTI ALFF*RSFA*; 3dAFNItoNIFTI ALFF*fRSFA*; 3dAFNItoNIFTI ALFF*mALFF*; 3dAFNItoNIFTI ALFF*mRSFA*
rm *.BRIK *.HEAD
gzip ALFF*

#Calculate ALFF stats
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -P 97.5)
ALFF4=$(fslstats ALFF_ALFF.nii.gz -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -P 97.5)
fALFF4=$(fslstats ALFF_fALFF.nii.gz -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -P 97.5)
mALFF4=$(fslstats ALFF_mALFF.nii.gz -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -P 97.5)
fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -P 97.5)
mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for sagsinus
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for CSF
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for GM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for WM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for allwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for dwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for pwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for nawm
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial only mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_ALFF4=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fALFF4=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mALFF4=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fRSFA4=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mRSFA4=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

cat >> ./"$ID"_alff4_data.txt << EOL
$1 $ALFF4 $fALFF4 $mALFF4 $fRSFA4 $mRSFA4 $arterial_ALFF4 $arterial_fALFF4 $arterial_mALFF4 $arterial_fRSFA4 $arterial_mRSFA4 $sagsinus_ALFF4 $sagsinus_fALFF4 $sagsinus_mALFF4 $sagsinus_fRSFA4 $sagsinus_mRSFA4 $CSF_ALFF4 $CSF_fALFF4 $CSF_mALFF4 $CSF_fRSFA4 $CSF_mRSFA4 $GM_ALFF4 $GM_fALFF4 $GM_mALFF4 $GM_fRSFA4 $GM_mRSFA4 $WM_ALFF4 $WM_fALFF4 $WM_mALFF4 $WM_fRSFA4 $WM_mRSFA4 $allwmh_ALFF4 $allwmh_fALFF4 $allwmh_mALFF4 $allwmh_fRSFA4 $allwmh_mRSFA4 $dwmh_ALFF4 $dwmh_fALFF4 $dwmh_mALFF4 $dwmh_fRSFA4 $dwmh_mRSFA4 $pwmh_ALFF4 $pwmh_fALFF4 $pwmh_mALFF4 $pwmh_fRSFA4 $pwmh_mRSFA4 $nawm_ALFF4 $nawm_fALFF4 $nawm_mALFF4 $nawm_fRSFA4 $nawm_mRSFA4 $arterial_only_ALFF4 $arterial_only_fALFF4 $arterial_only_mALFF4 $arterial_only_fRSFA4 $arterial_only_mRSFA4
EOL

cat ./"$ID"_alff4_data.txt >> "$homepath"final_outputs/all_alff4_data.txt

echo -------- $ID complete ALFF4 `date`

###################################################################################################################################################################################################

echo ----- $ID start ALFF0 `date`

cd "$homepath""$ID""_analysis/alff/alff0/"

##Create ALFF files
3dRSFC -prefix ALFF 0.01 0.05 "$homepath""$ID""_analysis/preprocessing/rfMRI_normalised.nii.gz"

##Organise ALFF files
3dAFNItoNIFTI ALFF*ALFF*; 3dAFNItoNIFTI ALFF*fALFF*; 3dAFNItoNIFTI ALFF*RSFA*; 3dAFNItoNIFTI ALFF*fRSFA*; 3dAFNItoNIFTI ALFF*mALFF*; 3dAFNItoNIFTI ALFF*mRSFA*
rm *.BRIK *.HEAD
gzip ALFF*

#Calculate ALFF stats
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -P 97.5)
ALFF0=$(fslstats ALFF_ALFF.nii.gz -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -P 97.5)
fALFF0=$(fslstats ALFF_fALFF.nii.gz -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -P 97.5)
mALFF0=$(fslstats ALFF_mALFF.nii.gz -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -P 97.5)
fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -P 97.5)
mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -P 97.5)
arterial_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_pulsatility_mask_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for sagsinus
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -P 97.5)
sagsinus_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_sagsinus_mod.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for CSF
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -P 97.5)
CSF_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_CSF_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for GM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -P 97.5)
GM_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_GM_exclusive.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for WM
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -P 97.5)
WM_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_brain_WM_exclusive_noDG.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for allwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -P 97.5)
allwmh_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_allwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for dwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -P 97.5)
dwmh_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_dwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for pwmh
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -P 97.5)
pwmh_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_bianca_pwmh.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for nawm
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -P 97.5)
nawm_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_nawm.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

#ALFF stats for arterial only mask
ALFF_robust_2p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
ALFF_robust_97p5=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_ALFF0=$(fslstats ALFF_ALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $ALFF_robust_2p5 -u $ALFF_robust_97p5 -M -S)

fALFF_robust_2p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fALFF_robust_97p5=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fALFF0=$(fslstats ALFF_fALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fALFF_robust_2p5 -u $fALFF_robust_97p5 -M -S)

mALFF_robust_2p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mALFF_robust_97p5=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mALFF0=$(fslstats ALFF_mALFF.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mALFF_robust_2p5 -u $mALFF_robust_97p5 -M -S)

fRSFA_robust_2p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
fRSFA_robust_97p5=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_fRSFA0=$(fslstats ALFF_fRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $fRSFA_robust_2p5 -u $fRSFA_robust_97p5 -M -S)

mRSFA_robust_2p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 2.5)
mRSFA_robust_97p5=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -P 97.5)
arterial_only_mRSFA0=$(fslstats ALFF_mRSFA.nii.gz -k "$homepath""$ID""_analysis/masks/BOLD_arterial_only.nii.gz" -M -S -P 50 -l $mRSFA_robust_2p5 -u $mRSFA_robust_97p5 -M -S)

cat >> ./"$ID"_alff0_data.txt << EOL
$1 $ALFF0 $fALFF0 $mALFF0 $fRSFA0 $mRSFA0 $arterial_ALFF0 $arterial_fALFF0 $arterial_mALFF0 $arterial_fRSFA0 $arterial_mRSFA0 $sagsinus_ALFF0 $sagsinus_fALFF0 $sagsinus_mALFF0 $sagsinus_fRSFA0 $sagsinus_mRSFA0 $CSF_ALFF0 $CSF_fALFF0 $CSF_mALFF0 $CSF_fRSFA0 $CSF_mRSFA0 $GM_ALFF0 $GM_fALFF0 $GM_mALFF0 $GM_fRSFA0 $GM_mRSFA0 $WM_ALFF0 $WM_fALFF0 $WM_mALFF0 $WM_fRSFA0 $WM_mRSFA0 $allwmh_ALFF0 $allwmh_fALFF0 $allwmh_mALFF0 $allwmh_fRSFA0 $allwmh_mRSFA0 $dwmh_ALFF0 $dwmh_fALFF0 $dwmh_mALFF0 $dwmh_fRSFA0 $dwmh_mRSFA0 $pwmh_ALFF0 $pwmh_fALFF0 $pwmh_mALFF0 $pwmh_fRSFA0 $pwmh_mRSFA0 $nawm_ALFF0 $nawm_fALFF0 $nawm_mALFF0 $nawm_fRSFA0 $nawm_mRSFA0 $arterial_only_ALFF0 $arterial_only_fALFF0 $arterial_only_mALFF0 $arterial_only_fRSFA0 $arterial_only_mRSFA0
EOL

cat ./"$ID"_alff0_data.txt >> "$homepath"final_outputs/all_alff0_data.txt

echo -------- $ID complete ALFF0 `date`


###################################################################################################################################################################################################

echo ----- $ID start monte carlo analysis `date`

cd "$homepath""$ID""_analysis/preprocessing/"

matlab -singleCompThread -nodisplay -nosplash -nojvm -r "run /well/webb/users/vox025/ViewPowerCalcs3($ID)"

cd "$homepath""$ID""_analysis/preprocessing/"

cat ./mc_analysis_values3.txt >> "$homepath"final_outputs/all_mc_analysis_values3.txt

cat ./transfer_analysis_values.txt >> "$homepath"final_outputs/all_transfer_analysis_values.txt

cat ./spectra_analysis_values.txt >> "$homepath"final_outputs/all_spectra_analysis_values.txt

#matlab -singleCompThread -nodisplay -nosplash -nojvm -r "run /well/webb/users/vox025/ViewPowerCalcs0($ID)"

#cd "$homepath""$ID""_analysis/preprocessing/"

#cat ./mc_analysis_values0.txt >> "$homepath"final_outputs/all_mc_analysis_values0.txt

matlab -singleCompThread -nodisplay -nosplash -nojvm -r "run /well/webb/users/vox025/ViewPowerCalcs1($ID)"

cd "$homepath""$ID""_analysis/preprocessing/"

cat ./mc_analysis_values1.txt >> "$homepath"final_outputs/all_mc_analysis_values1.txt

matlab -singleCompThread -nodisplay -nosplash -nojvm -r "run /well/webb/users/vox025/ViewPowerCalcs2($ID)"

cd "$homepath""$ID""_analysis/preprocessing/"

cat ./mc_analysis_values2.txt >> "$homepath"final_outputs/all_mc_analysis_values2.txt

#matlab -singleCompThread -nodisplay -nosplash -nojvm -r "run /well/webb/users/vox025/ViewPowerCalcs4($ID)"

#cd "$homepath""$ID""_analysis/preprocessing/"

#cat ./mc_analysis_values4.txt >> "$homepath"final_outputs/all_mc_analysis_values4.txt

echo -------- $ID complete monte carlo analysis `date`

###################################################################################################################################################################################################

#echo ----- $ID start final data file `date`

#cd "$homepath""$ID""_analysis/preprocessing/"

#cat >> ./"$ID"_all_data.txt << EOL
#$1 $T1 $T2_FLAIR $rfMRI $SWI $rfMRI_pspec_in_standard $ALFF0_in_standard $ALFF1_in_standard $ALFF2_in_standard $ALFF3_in_standard $ALFF4_in_standard $masked_ts_arterial $masked_ts_sagsinus $masked_ts_CSF $masked_ts_GM $masked_ts_WM $masked_ts_allwmh $masked_ts_dwmh $masked_ts_pwmh $masked_ts_nawm $masked_ts_arterial_only $global_pspec $arterial_pspec $sagsinus_pspec $CSF_pspec $GM_pspec $WM_pspec $allwmh_pspec $dwmh_pspec $pwmh_pspec $nawm_pspec $arterial_only_pspec $ALFF0 $fALFF0 $mALFF0 $fRSFA0 $mRSFA0 $arterial_ALFF0 $arterial_fALFF0 $arterial_mALFF0 $arterial_fRSFA0 $arterial_mRSFA0 $sagsinus_ALFF0 $sagsinus_fALFF0 $sagsinus_mALFF0 $sagsinus_fRSFA0 $sagsinus_mRSFA0 $CSF_ALFF0 $CSF_fALFF0 $CSF_mALFF0 $CSF_fRSFA0 $CSF_mRSFA0 $GM_ALFF0 $GM_fALFF0 $GM_mALFF0 $GM_fRSFA0 $GM_mRSFA0 $WM_ALFF0 $WM_fALFF0 $WM_mALFF0 $WM_fRSFA0 $WM_mRSFA0 $allwmh_ALFF0 $allwmh_fALFF0 $allwmh_mALFF0 $allwmh_fRSFA0 $allwmh_mRSFA0 $dwmh_ALFF0 $dwmh_fALFF0 $dwmh_mALFF0 $dwmh_fRSFA0 $dwmh_mRSFA0 $pwmh_ALFF0 $pwmh_fALFF0 $pwmh_mALFF0 $pwmh_fRSFA0 $pwmh_mRSFA0 $nawm_ALFF0 $nawm_fALFF0 $nawm_mALFF0 $nawm_fRSFA0 $nawm_mRSFA0 $arterial_only_ALFF0 $arterial_only_fALFF0 $arterial_only_mALFF0 $arterial_only_fRSFA0 $arterial_only_mRSFA0 $ALFF1 $fALFF1 $mALFF1 $fRSFA1 $mRSFA1 $arterial_ALFF1 $arterial_fALFF1 $arterial_mALFF1 $arterial_fRSFA1 $arterial_mRSFA1 $sagsinus_ALFF1 $sagsinus_fALFF1 $sagsinus_mALFF1 $sagsinus_fRSFA1 $sagsinus_mRSFA1 $CSF_ALFF1 $CSF_fALFF1 $CSF_mALFF1 $CSF_fRSFA1 $CSF_mRSFA1 $GM_ALFF1 $GM_fALFF1 $GM_mALFF1 $GM_fRSFA1 $GM_mRSFA1 $WM_ALFF1 $WM_fALFF1 $WM_mALFF1 $WM_fRSFA1 $WM_mRSFA1 $allwmh_ALFF1 $allwmh_fALFF1 $allwmh_mALFF1 $allwmh_fRSFA1 $allwmh_mRSFA1 $dwmh_ALFF1 $dwmh_fALFF1 $dwmh_mALFF1 $dwmh_fRSFA1 $dwmh_mRSFA1 $pwmh_ALFF1 $pwmh_fALFF1 $pwmh_mALFF1 $pwmh_fRSFA1 $pwmh_mRSFA1 $nawm_ALFF1 $nawm_fALFF1 $nawm_mALFF1 $nawm_fRSFA1 $nawm_mRSFA1 $arterial_only_ALFF1 $arterial_only_fALFF1 $arterial_only_mALFF1 $arterial_only_fRSFA1 $arterial_only_mRSFA1 $ALFF2 $fALFF2 $mALFF2 $fRSFA2 $mRSFA2 $arterial_ALFF2 $arterial_fALFF2 $arterial_mALFF2 $arterial_fRSFA2 $arterial_mRSFA2 $sagsinus_ALFF2 $sagsinus_fALFF2 $sagsinus_mALFF2 $sagsinus_fRSFA2 $sagsinus_mRSFA2 $CSF_ALFF2 $CSF_fALFF2 $CSF_mALFF2 $CSF_fRSFA2 $CSF_mRSFA2 $GM_ALFF2 $GM_fALFF2 $GM_mALFF2 $GM_fRSFA2 $GM_mRSFA2 $WM_ALFF2 $WM_fALFF2 $WM_mALFF2 $WM_fRSFA2 $WM_mRSFA2 $allwmh_ALFF2 $allwmh_fALFF2 $allwmh_mALFF2 $allwmh_fRSFA2 $allwmh_mRSFA2 $dwmh_ALFF2 $dwmh_fALFF2 $dwmh_mALFF2 $dwmh_fRSFA2 $dwmh_mRSFA2 $pwmh_ALFF2 $pwmh_fALFF2 $pwmh_mALFF2 $pwmh_fRSFA2 $pwmh_mRSFA2 $nawm_ALFF2 $nawm_fALFF2 $nawm_mALFF2 $nawm_fRSFA2 $nawm_mRSFA2 $arterial_only_ALFF2 $arterial_only_fALFF2 $arterial_only_mALFF2 $arterial_only_fRSFA2 $arterial_only_mRSFA2 $ALFF3 $fALFF3 $mALFF3 $fRSFA3 $mRSFA3 $arterial_ALFF3 $arterial_fALFF3 $arterial_mALFF3 $arterial_fRSFA3 $arterial_mRSFA3 $sagsinus_ALFF3 $sagsinus_fALFF3 $sagsinus_mALFF3 $sagsinus_fRSFA3 $sagsinus_mRSFA3 $CSF_ALFF3 $CSF_fALFF3 $CSF_mALFF3 $CSF_fRSFA3 $CSF_mRSFA3 $GM_ALFF3 $GM_fALFF3 $GM_mALFF3 $GM_fRSFA3 $GM_mRSFA3 $WM_ALFF3 $WM_fALFF3 $WM_mALFF3 $WM_fRSFA3 $WM_mRSFA3 $allwmh_ALFF3 $allwmh_fALFF3 $allwmh_mALFF3 $allwmh_fRSFA3 $allwmh_mRSFA3 $dwmh_ALFF3 $dwmh_fALFF3 $dwmh_mALFF3 $dwmh_fRSFA3 $dwmh_mRSFA3 $pwmh_ALFF3 $pwmh_fALFF3 $pwmh_mALFF3 $pwmh_fRSFA3 $pwmh_mRSFA3 $nawm_ALFF3 $nawm_fALFF3 $nawm_mALFF3 $nawm_fRSFA3 $nawm_mRSFA3 $arterial_only_ALFF3 $arterial_only_fALFF3 $arterial_only_mALFF3 $arterial_only_fRSFA3 $arterial_only_mRSFA3 $ALFF4 $fALFF4 $mALFF4 $fRSFA4 $mRSFA4 $arterial_ALFF4 $arterial_fALFF4 $arterial_mALFF4 $arterial_fRSFA4 $arterial_mRSFA4 $sagsinus_ALFF4 $sagsinus_fALFF4 $sagsinus_mALFF4 $sagsinus_fRSFA4 $sagsinus_mRSFA4 $CSF_ALFF4 $CSF_fALFF4 $CSF_mALFF4 $CSF_fRSFA4 $CSF_mRSFA4 $GM_ALFF4 $GM_fALFF4 $GM_mALFF4 $GM_fRSFA4 $GM_mRSFA4 $WM_ALFF4 $WM_fALFF4 $WM_mALFF4 $WM_fRSFA4 $WM_mRSFA4 $allwmh_ALFF4 $allwmh_fALFF4 $allwmh_mALFF4 $allwmh_fRSFA4 $allwmh_mRSFA4 $dwmh_ALFF4 $dwmh_fALFF4 $dwmh_mALFF4 $dwmh_fRSFA4 $dwmh_mRSFA4 $pwmh_ALFF4 $pwmh_fALFF4 $pwmh_mALFF4 $pwmh_fRSFA4 $pwmh_mRSFA4 $nawm_ALFF4 $nawm_fALFF4 $nawm_mALFF4 $nawm_fRSFA4 $nawm_mRSFA4 $arterial_only_ALFF4 $arterial_only_fALFF4 $arterial_only_mALFF4 $arterial_only_fRSFA4 $arterial_only_mRSFA4
#EOL

#cat ./"$ID"_all_data.txt >> "$homepath"final_outputs/all_participants_all_data.txt

#echo -------- $ID complete final data file `date`

###################################################################################################################################################################################################

cd "$homepath"
rm -r "$ID""_analysis"/

echo -------- $ID complete all processing `date`
