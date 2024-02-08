#!/usr/bin/env Rscript

#edited by Longchuan--------
#library("refund",lib.loc="/home/lli/R/x86_64-redhat-linux-gnu-library/3.2/")
library("amap")
library("mgcv")
require("ggplot2")
library("data.table")
library("itsadug")
library("tools")
library("gridExtra")
require("MASS")
library("pspline")
library("reshape2")
library("dendextend")
library("ggdendro")
library("graphics")
library("raster")

#dev.off()

path="/home/lli/projects_dti_harmonization_for_site_by_BRisk/Pipelines_ExampleData/hierarch_fatr_fine_1mm_mode3/init_temp/final_temp/diffeo_warped"
#filename_input="dat_harmonized_from_BRisk_WB_old_format.txt"
filename_input="/dat_harmonized_from_BRisk_11Tracts_old_format.txt"

#path="/home/lli/projects_TnP2020_infant_white_matter_development_updated_09_08_2020/Pipelines_ExampleData/hierarch_fatr_fine_1mm_mode3/init_temp/final_temp/diffeo_warped"
#filename_input="/dti_data_extracted_by_tracts_all_atlas_symmetry_xyflipped_2_TnP2020.txt"
filename=paste(path, filename_input, sep="")
base_name=file_path_sans_ext(basename(filename_input));
filename_opt=paste(path,"/",base_name,"_gamm.csv", sep="")
if (file.exists(filename_opt)){
file.remove(filename_opt)
}

#deri_method='Diff' # Diff: use diff() to derive derivatives; Pspline: using Smooth.Pspline function
fit_method="REML"
mah_type=1 # mahalanobis distance type: 1: FA+RD; 2: FA+RD+AD; 3: RD+AD; 4: FA+MD+AD+RD

bs_val='tp' # tp, cr, bs, 
norder_val=4 # number of order in smoothing the curves
df_val=12
m_val=3
k_val=20
gamma_val=0.5
sp_val=c(1)
select_pen=FALSE
nrep=20
nt=200 # number of time points 
cr_type=2 #1: original change rate, 2: percentile change rate

#filename_optS=paste(path,"/",base_name,"S.csv", sep="")
dir_name=dirname(filename);
mac=read.table(filename, header=TRUE)
uniDTI=unique(mac$DTI_type)
uniReg=unique(mac$ROI)
uniReg_len=length(uniReg)

hcp_filename="/home/lli/projects_adults_ref_mac_hcp_combined/Pipelines_ExampleData/hierarch_fatr_fine_1mm_mode3/init_temp/final_temp/diffeo_warped/dti_data_extracted_by_tracts_all_atlas_symmetry_xyflipped_2_fa_mean_adults.txt"
hcp=read.table(hcp_filename, header=TRUE)

col_tab=read.table(file="/home/lli/matlab/Myscript_Unix/Myscripts_shared/table_included_tracts_reorderedCol_sym_hex.txt", header=FALSE, colClasses="character")

# original model terms
smth_gam_all=array(list(), c(length(uniDTI), length(uniReg)))
interp_gam_all=smth_gam_all

temp=mac[mac$DTI_type==uniDTI[1] & mac$Status=="LR" & mac$ROI==uniReg[1],];
mac_mtx=array(NaN,c(length(uniDTI), length(uniReg), length(temp$DTI_val)))

temp=hcp[hcp$DTI_type==uniDTI[1] & hcp$ROI==uniReg[1],];
hcp_mtx=array(NaN,c(length(uniDTI), length(uniReg), length(temp$DTI_val)))

for (i in 1:length(uniDTI))
{
	for (j in 1:length(uniReg)) 
	{
		print(paste0(uniDTI[i],' : ',uniReg[j]));
		rm(macT, hcpT)
		macT=mac[mac$DTI_type==uniDTI[i] & mac$Status=="LR" & mac$ROI==uniReg[j],];
		hcpT=hcp[hcp$DTI_type==uniDTI[i] & hcp$ROI==uniReg[j],]
		mac_mtx[i,j,]=c(macT$DTI_val)
		hcp_mtx[i,j,]=c(hcpT$DTI_val)
	}
}

color_tab=rep(paste0("#",col_tab[,8]), each=length(macT$Corr_age))

print("calculating Mahalanobis distance ...")
# mahalanobis distance type: 1: FA+RD; 2: FA+RD+AD; 3: RD+AD; 4: FA+MD+AD+RD; 5: FA; 6: MD; 7: AD; 8: RD
mac_dim=dim(mac_mtx)
df_comb=data.frame()
for (mah_type in 1:3)
{
	rm(df_all)
	df_all=data.frame()
	for (j in 1:length(uniReg))
	{
		rm(x,center,mah_dist, x, y)
		mah_dist=array(NaN,c(3, mac_dim[2], mac_dim[3]))
		if ( mah_type == 1)
		{
			#FA+RD
			x=cbind(mac_mtx[1,j,],mac_mtx[4,j,])
			y=cbind(hcp_mtx[1,j,], hcp_mtx[4,j,])
			
		} else if ( mah_type ==2)
		{
			# FA+RD+AD
			x=cbind(mac_mtx[1,j,],mac_mtx[3,j,], mac_mtx[4,j,])
			y=cbind(hcp_mtx[1,j,], hcp_mtx[3,j,], hcp_mtx[4,j,])
			
		} else if (mah_type == 3)
		{
			# RD+AD
			x=cbind(mac_mtx[3,j,], mac_mtx[4,j,])
			y=cbind(hcp_mtx[3,j,], hcp_mtx[4,j,])
		} else if (mah_type ==4) 
		{
			x=cbind(mac_mtx[1,j,], mac_mtx[2,j,],mac_mtx[3,j,], mac_mtx[4,j,])
			y=cbind(hcp_mtx[1,j,], hcp_mtx[2,j,], hcp_mtx[3,j,], hcp_mtx[4,j,])
			
		} else if (mah_type ==5)
		{
			x=cbind(mac_mtx[1,j,])
			y=cbind(hcp_mtx[1,j,])
		}
		
		 else if (mah_type ==6)
		{
			x=cbind(mac_mtx[2,j,])
			y=cbind(hcp_mtx[2,j,])
		}
		
		 else if (mah_type ==7)
		{
			x=cbind(mac_mtx[3,j,])
			y=cbind(hcp_mtx[3,j,])
		}
		
		 else if (mah_type ==8)
		{
			x=cbind(mac_mtx[4,j,])
			y=cbind(hcp_mtx[4,j,])
		}
		
		center=colMeans(y)
		mah_dist[mah_type,j,]=mahalanobis(x, center, var(y))
		
		macT=mac[mac$DTI_type==uniDTI[1] & mac$Status=="LR" & mac$ROI==uniReg[j],];
		temp_df <- data.frame(Unique_ID=macT$Unique_ID, mah_type=rep(mah_type, each=mac_dim[3]), Corr_age=macT$Corr_age, mah_dis=mah_dist[mah_type, j,], color_tab=rep(paste0("#",col_tab[j,8]), each=mac_dim[3]), region=as.factor(rep(uniReg[j], each=mac_dim[3])))
		df_all <- rbind(df_all,temp_df)
	}
	df_comb=rbind(df_comb, df_all)
	
	x_label=paste("Age (days)")
	xlim_val=c(-10, 200)
	y_label="Mah. Dist."
}


print("fit the mahalanobis distance curve...")
# original model terms
smth_gam_all=array(list(), c(length(uniDTI), length(uniReg)))
interp_gam_all=smth_gam_all
df_zero=data.frame()

# fit the curve of the mahalanobis distance
df_mah_comb=data.frame()
for (mah_type in 1:length(unique(df_comb$mah_type)))
{
	rm(df_all)
	df_all=data.frame()
	rm(zero_t)
	zero_t=array(data=NA, dim=length(uniReg))
	for (j in 1:mac_dim[2])
	{
		mahT=df_comb[df_comb$mah_type==mah_type & df_comb$region==uniReg[j],];
		x=seq(from=min(mahT$Corr_age), to=max(mahT$Corr_age), length.out=nt)
		x_ext=x
		paste0("x:",x)
		paste0("x_ext:",x_ext)

		#m5_gam=gam(mah_dis~ s(Corr_age, k=k_val, m=m_val, bs=bs_val) , method=fit_method, data=mahT, gamma=gamma_val, select=select_pen, sp=sp_val)
		
		m5_gam=gam(mah_dis~ s(Corr_age)+s(Unique_ID,bs='re'),method='REML',data=dat_harmonized)
		
		smth_gam <- get_modelterm(m5_gam, select=1, n.grid=nt, se=TRUE)
		interp_gam <- get_coefs(m5_gam)
		smth_gam_all[[mah_type,j]] <- smth_gam
		interp_gam_all[[mah_type,j]] <- get_coefs(m5_gam)
		x_label=pastte("Age (days)")
		xlim_val=c(-10, 200)
		y_label="Mah. Dist."
		
		y_int=predict(sm.spline(smth_gam_all[[mah_type,j]]$terms, smth_gam_all[[mah_type,j]]$fit+interp_gam_all[[mah_type,j]][[1]], norder=norder_val, df=df_val), x_ext)		
		zero_orig_x=array(data=NA, dim=length(x))
		zero_orig_y=zero_orig_x

		temp_df <- data.frame(mah_type=mah_type, age_corr_new=x, y_val=smth_gam$fit + interp_gam[[1]], color_tab=rep(paste0("#",col_tab[j,8]), each=nt), region=as.factor(rep(uniReg[j], each=nt)))
		df_all <- rbind(df_all,temp_df)
	}
	df_mah_comb=rbind(df_mah_comb, df_all)
	df_tmp=data.frame(t(zero_t))
	df_zero<-rbind(df_zero, df_tmp)
}

print("fitting the first derivative (dr1%) of mahalanobis distance...")
# plot 1st derivatives (%)
df_mah_per_comb=data.frame()
df_mah_per_orig_comb=data.frame()

for (mah_type in 1:length(unique(df_comb$mah_type)))
{
	rm(df_all, zero_t)
	df_all=data.frame()
	df_all_orig=data.frame()
	zero_t=array(data=NA, dim=length(uniReg))
	for (j in 1:length(uniReg)) 
	{
		print(paste0(' : ',uniReg[j]));

		tmp0<-predict(sm.spline(smth_gam_all[[mah_type,j]]$terms, smth_gam_all[[mah_type,j]]$fit+interp_gam_all[[mah_type,j]][[1]], norder=norder_val, df=df_val), x_ext)
		tmp1<-predict(sm.spline(smth_gam_all[[mah_type,j]]$terms, smth_gam_all[[mah_type,j]]$fit+interp_gam_all[[mah_type,j]][[1]], norder=norder_val, df=df_val), x_ext,1)
		y_int_dr1_per=100*tmp1/tmp0
		y_int_dr1_per_orig=tmp1

		x_label=paste("Age (days)")
		xlim_val=c(-10, 200)
		y_label="Mah. Dist. C.R. (%)"

		#generating a data frame for plotting
		temp_df <- data.frame(mah_type=mah_type, age_corr_new=x, y_val=y_int_dr1_per, color_tab=rep(paste0("#",col_tab[j,8]), each=nt), region=as.factor(rep(uniReg[j], each=nt)))
		df_all <- rbind(df_all,temp_df)
		
		temp_df_orig <- data.frame(mah_type=mah_type, age_corr_new=x, y_val=y_int_dr1_per_orig, color_tab=rep(paste0("#",col_tab[j,8]), each=nt), region=as.factor(rep(uniReg[j], each=nt)))
		df_all_orig <- rbind(df_all_orig,temp_df_orig)
	}

	df_mah_per_comb=rbind(df_mah_per_comb, df_all)
	df_mah_per_orig_comb=rbind(df_mah_per_orig_comb, df_all_orig)
	df_tmp=data.frame(t(zero_t))
	df_zero<-rbind(df_zero, df_tmp)
}

	
print("identify the relationships between absolute mahalanobis distance and its dynamics (%/day)...")
age_corr_new_uni=unique(df_mah_comb$age_corr_new)
stats_test=data.frame()
for (mah_type in 1:length(unique(df_comb$mah_type)))
{
	rm(df_t)
	df_t=data.frame()
	for (t in 1: length(age_corr_new_uni))
	{
		rm(tmp_orig,tmp_per)
		tmp_orig=df_mah_comb[df_mah_comb$mah_type==mah_type & df_mah_comb$age_corr_new==age_corr_new_uni[t],]
		tmp_per=df_mah_per_comb[df_mah_per_comb$mah_type==mah_type & df_mah_per_comb$age_corr_new==age_corr_new_uni[t],]
        rm(s_tmp, time, pearson_r, pearson_cidn, pearson_ciup, spearman_r, mah, region)
        time = age_corr_new_uni[t]
        pearson_r=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("pearson"))$estimate
        pearson_cidn=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("pearson"))$conf.int[1]
        pearson_ciup=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("pearson"))$conf.int[2]
        spearman_r=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("spearman"))$estimate
        pearson_p=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("pearson"))$p.value
        mah=mah_type
        region=uniReg[j]
        s_tmp<-data.frame(time=time, pearson_p=pearson_p, pearson_r=pearson_r, pearson_cidn=pearson_cidn, pearson_ciup=pearson_ciup, spearman_r=spearman_r, mah=mah, region=region)
        df_t=rbind(df_t, s_tmp)
	}
	stats_test=rbind(stats_test,df_t)
	nam<-paste("p",mah_type, sep="")
	coln<- y_label
	attach(df_all)
	assign(
	nam, ggplot(data=df_t) +
	geom_line( aes_(x=df_t$time, y=df_t$pearson_r, color="black"),  size=2, show.legend=FALSE) +
	geom_ribbon(aes_(x=df_t$time, ymin=df_t$pearson_cidn, ymax=df_t$pearson_ciup), fill="gray", alpha=0.5) +
	#geom_line(aes_(x=df_t$time, y=df_t$pearson_p, color="red"), size=2, show.legend=FALSE) +
	labs(x=x_label, y="Correlation Coefficient") +
	theme(axis.text=element_text(size=16), axis.title = element_text(size = 18)) 
	+geom_hline(yintercept=0, color='black', size=1)
	)
}

pdf(file=paste(path,"/",base_name,"_mah_corr_bet_orig_and_speed_all.pdf", sep=""), width=15, height=4)
plist<-mget(paste0("p", 1:3))
do.call(grid.arrange, c(plist,nrow=1))
dev.off()



print("identify the relationships between absolute mahalanobis distance and its dynamics (original values /per day...")
age_corr_new_uni=unique(df_mah_comb$age_corr_new)
stats_test=data.frame()
for (mah_type in 1:length(unique(df_comb$mah_type)))
{
	rm(df_t)
	df_t=data.frame()
	for (t in 1: length(age_corr_new_uni))
	{
		rm(tmp_orig,tmp_per)
		tmp_orig=df_mah_comb[df_mah_comb$mah_type==mah_type & df_mah_comb$age_corr_new==age_corr_new_uni[t],]
		tmp_per=df_mah_per_orig_comb[df_mah_per_orig_comb$mah_type==mah_type & df_mah_per_orig_comb$age_corr_new==age_corr_new_uni[t],]
        rm(s_tmp, time, pearson_r, pearson_cidn, pearson_ciup, spearman_r, mah, region)
        time = age_corr_new_uni[t]
        pearson_r=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("pearson"))$estimate
        pearson_cidn=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("pearson"))$conf.int[1]
        pearson_ciup=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("pearson"))$conf.int[2]
        spearman_r=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("spearman"))$estimate
        pearson_p=cor.test(tmp_orig$y_val, tmp_per$y_val, method=c("pearson"))$p.value
        mah=mah_type
        region=uniReg[j]
        s_tmp<-data.frame(time=time, pearson_p=pearson_p, pearson_r=pearson_r, pearson_cidn=pearson_cidn, pearson_ciup=pearson_ciup, spearman_r=spearman_r, mah=mah, region=region)
        df_t=rbind(df_t, s_tmp)
	}
	stats_test=rbind(stats_test,df_t)
	nam<-paste("p",mah_type, sep="")
	coln<- y_label
	attach(df_all)
	assign(
	nam, ggplot(data=df_t) +
	geom_line( aes_(x=df_t$time, y=df_t$pearson_r, color="black"),  size=2, show.legend=FALSE) +
	geom_ribbon(aes_(x=df_t$time, ymin=df_t$pearson_cidn, ymax=df_t$pearson_ciup), fill="gray", alpha=0.5) +
	#geom_line(aes_(x=df_t$time, y=df_t$pearson_p, color="red"), size=2, show.legend=FALSE) +
	labs(x=x_label, y="Correlation Coefficient") +
	theme(axis.text=element_text(size=16), axis.title = element_text(size = 18)) 
	+geom_hline(yintercept=0, color='black', size=1)
	)
}

pdf(file=paste(path,"/",base_name,"_mah_corr_bet_orig_and_speed_orig_all.pdf", sep=""), width=15, height=4)
plist<-mget(paste0("p", 1:3))
do.call(grid.arrange, c(plist,nrow=1))
dev.off()


print(paste0("leaving mac_infants_humans_mahalanobis_distance.R...", sep=""))

