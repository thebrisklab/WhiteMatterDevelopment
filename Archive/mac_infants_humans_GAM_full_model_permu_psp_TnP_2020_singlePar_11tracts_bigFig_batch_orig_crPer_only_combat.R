#!/usr/bin/env Rscript

#edited by Longchuan--------
#library("refund",lib.loc="/home/lli/R/x86_64-redhat-linux-gnu-library/3.2/")

library("mgcv")
require("ggplot2")
library("data.table")
library("itsadug")
library("tools")
library("gridExtra")
require("MASS")
library("pspline")
library("reshape2")
#dev.off()

path="/home/lli/projects_dti_harmonization_for_site_by_BRisk/Pipelines_ExampleData/hierarch_fatr_fine_1mm_mode3/init_temp/final_temp/diffeo_warped"
#filename_input="dat_harmonized_from_BRisk_WB_old_format.txt"
filename_input="/dat_harmonized_from_BRisk_11Tracts_old_format.txt"


filename=paste(path,"/", filename_input, sep="")
base_name=file_path_sans_ext(basename(filename_input));
filename_opt=paste(path,"/",base_name,"_gamm.csv", sep="")
if (file.exists(filename_opt)){
file.remove(filename_opt)
}

#deri_method='Diff' # Diff: use diff() to derive derivatives; Pspline: using Smooth.Pspline function
ifRelative=1 # 0: original plot; 1: relative white matter devleopment normalized by the mean; 2: relative white matter development relative how many SD  
cr_per=1 # 0: plot original change rate curve, 1: plot percentile change rates, which is easier to interpret

if (ifRelative > 0)
{
	hcp_filename="/home/lli/projects_adults_ref_mac_hcp_combined/Pipelines_ExampleData/hierarch_fatr_fine_1mm_mode3/init_temp/final_temp/diffeo_warped/dti_data_extracted_by_tracts_all_atlas_symmetry_xyflipped_2_fa_mean_adults.txt"
	hcp=read.table(hcp_filename, header=TRUE)
}

fit_method="REML"
#sp_val=c(20)
bs_val='tp' # tp, cr, bs, 
norder_val=4 # number of order in smoothing the curves
norder_val=3 # testing new values
m_val=3
k_val=20
gamma_val=0.5
df_val=12
sp_val=c(1)
select_pen=FALSE
nrep=500
nt=200 # number of time points 
#filename_optS=paste(path,"/",base_name,"S.csv", sep="")
dir_name=dirname(filename);
mac=read.table(filename, header=TRUE)
uniDTI=unique(mac$DTI_type)
uniReg=unique(mac$ROI)
uniReg_len=length(uniReg)

col_tab=read.table(file="/home/lli/matlab/Myscript_Unix/Myscripts_shared/table_included_tracts_reorderedCol_sym_hex.txt", header=FALSE, colClasses="character")

smth.b.all_sp_org=array(numeric(),dim = c(length(uniDTI), length(uniReg), nt, nrep))
smth.b.all_sp_dr1=smth.b.all_sp_org
smth.b.all_sp_dr2=smth.b.all_sp_org

# original model terms
smth_gam_all=array(list(), c(length(uniDTI), length(uniReg)))
interp_gam_all=smth_gam_all

filename_opt_csv=paste(path,"/",base_name,"_tracts_to_zero.csv", sep="")
if (file.exists(filename_opt_csv) == TRUE)
{file.remove(filename_opt_csv) }
df_zero=data.frame()

df_comb=data.frame() # a data frame to contain all the fitted values
for (i in 1:length(uniDTI))
{
	rm(df_all)
	df_all=data.frame()
	rm(zero_t)
	zero_t=array(data=NA, dim=length(uniReg))
	for (j in 1:length(uniReg)) 
	{
		print(paste0(uniDTI[i],' : ',uniReg[j]));
		rm(macT)
		macT=mac[mac$DTI_type==uniDTI[i] & mac$Status=="LR" & mac$ROI==uniReg[j],];
		macT$OFScanner=as.factor(ifelse(grepl('CSIP', macT$Actual_ID), 'Prisma', 'Trio'));
		macT$Scanner=ifelse(grepl('CSIP', macT$Actual_ID), 'Prisma', 'Trio');
		macT$OFScanner = as.ordered(macT$OFScanner);
		contrasts(macT$OFScanner) = 'contr.treatment'
		macT$Subject=as.factor(macT$Unique_ID);
		macT$OFSex = as.ordered(macT$Sex);
		contrasts(macT$OFSex) = 'contr.treatment'
		macT$OFGA = as.ordered(macT$GA);
		contrasts(macT$OFGA) = 'contr.treatment'
		
		hcpT=hcp[hcp$DTI_type==uniDTI[i] & hcp$ROI==uniReg[j],]
			
				# convert TR to MD
		if ( as.character(uniDTI[i]) == 'tr')
		{ 
			macT$DTI_val = macT$DTI_val/3
			hcpT$DTI_val = hcpT$DTI_val/3
			y_label="MD" 
		}
		
		# check whether adult reference need to be used to normalize
		if (ifRelative==1)
		{
			macT$DTI_val = macT$DTI_val/mean(hcpT$DTI_val)
		} else if (ifRelative==2)
		{
			macT$DTI_val = abs(macT$DTI_val  - mean(hcpT$DTI_val))/sd(hcpT$DTI_val)
	    }
		
		# reshuffle the data for construct bootstrapped CI at 0.05
		smth.b.all=array(list(), nrep)
		interp.b.all=array(list(), nrep)
		for (n in 1:nrep)
		{
			indiv <- unique(macT$Unique_ID)
			smp <- as.character(sort(sample(indiv, length(indiv), replace=TRUE)))
			smp.df <- data.frame(id=smp)
			macT.b = merge(smp.df, macT, by.x="id", by.y="Unique_ID")    
			n.obs <- table(macT$Unique_ID)
			smpU <- unique(smp)
			reps <- as.vector(table(smp))
			obs <- as.vector(n.obs[match(smpU, names(n.obs))])
			id.rep.obs <- cbind(as.numeric(gsub("_",".",smpU)), as.integer(reps), as.integer(obs))   
			NameFun <- function(info) {
			names <- paste0(rep(info[1], info[2]), ".", seq(1, info[2]))
			names.long <- sort(rep(names, info[3]))
			}

			macT.b[, 1] <- do.call("c", apply(id.rep.obs, 1, NameFun))
			#macT.b <- macT.b[order(macT.b[, 1], macT.b$Corr_age), ]
			macT.b$Unique_ID = macT.b$id
			#print(paste0("Unique ID #: ", length(unique(macT.b$Unique_ID))))
			m5.b=gam(DTI_val~ s(Corr_age, k=k_val, m=m_val, bs=bs_val) , method=fit_method, data=macT.b, gamma=gamma_val, select=select_pen, sp=sp_val)
			smth.b.all[[n]] <- get_modelterm(m5.b, select=1, n.grid=nt, se=FALSE, print.summary=FALSE)
			interp.b.all[[n]]<-get_coefs(m5.b)
		}

		x=seq(from=min(macT$Corr_age), to=max(macT$Corr_age), length.out=nt)
		x_knots=seq(from=min(x), to=max(x), length.out=k_val)

		#x_ext=seq(from=min(x_knots), to=max(x_knots), length.out=nt)
		x_ext=x
		paste0("x:",x)
		paste0("x_knots:",x_knots)
		paste0("x_ext:",x_ext)

		m5_gam=gam(DTI_val~ s(Corr_age, k=k_val, m=m_val, bs=bs_val) , knots=list(x_knots), method=fit_method, data=macT, gamma=gamma_val, select=select_pen, sp=sp_val)
		smth_gam <- get_modelterm(m5_gam, select=1, n.grid=nt, se=TRUE)
		interp_gam <- get_coefs(m5_gam)
		smth_gam_all[[i,j]] <- smth_gam
		interp_gam_all[[i,j]] <- get_coefs(m5_gam)

		#plot(smth_gam$terms, smth_gam$fit+interp_gam[1],"l", col="red")
		for (n in 1:nrep){
			#lines(smth.b.all[[n]]$terms, smth.b.all[[n]]$fit+interp.b.all[[n]][[1]],"l")
			tmp0=predict(smooth.Pspline(smth.b.all[[n]]$terms, smth.b.all[[n]]$fit+interp.b.all[[n]][[1]], norder=norder_val, df=df_val), x)
			smth.b.all_sp=tmp0
			# if ( anyNA(tmp1) == TRUE)  { print(paste0("nan found for n=",n," bootstrapping!")) }
			smth.b.all_sp_org[i,j,,n]=matrix(smth.b.all_sp)
			
			tmp1=predict(smooth.Pspline(smth.b.all[[n]]$terms, smth.b.all[[n]]$fit+interp.b.all[[n]][[1]], norder=norder_val, df=df_val), x,1)
			
			if (cr_per == 0)
			{ smth.b.all_sp_dr1[i,j,,n]=matrix(tmp1)
			} else {
			smth.b.all_sp_dr1[i,j,,n]=matrix(100*tmp1/tmp0)
			}

			tmp2=predict(smooth.Pspline(smth.b.all[[n]]$terms, smth.b.all[[n]]$fit+interp.b.all[[n]][[1]], norder=norder_val, df=df_val), x,2)
			smth.b.all_sp_dr2[i,j,,n]=matrix(tmp2)
			rm(tmp0, tmp1, tmp2)
		}

		op=smth.b.all_sp_org[i,j,,]
		tmp <- apply(op, 1, quantile, na.rm=TRUE, probs = c(0.025, 0.975))
		orig.ci=(tmp[2,] - tmp[1,])/2

		x_label=paste("Age (days)")
		xlim_val=c(-10, 200)
		y_label=toupper(uniDTI[i])

		if (as.character(uniDTI[i]) == "tr"){ y_label="MD" }
		if ( ifRelative ==0)
		{
			if ( i == 1) {ylim_val=c(0.1, 0.3) }  else {  ylim_val=c(0.75, 1.75) }
		} else  if ( ifRelative ==1 )
		{
			if ( i == 1) {ylim_val=c(0.25, 0.8) }  else {  ylim_val=c(0.9, 2.1) }
		}  else if ( ifRelative ==2 ) 
		{
			ylim_val=c(-0, 40) 
		}
	
		y_int=predict(sm.spline(smth_gam_all[[i,j]]$terms, smth_gam_all[[i,j]]$fit+interp_gam_all[[i,j]][[1]], norder=norder_val, df=df_val), x_ext)		
		zero_orig_x=array(data=NA, dim=length(x))
		zero_orig_y=zero_orig_x
		ciup=y_int + orig.ci
		cidn=y_int - orig.ci
		for ( t in 1:length(x))
		{
			if ((ciup[t]>0) & (cidn[t]<0))
			{
				zero_t[j]=x[t]
				break
			}	
		}

		print(paste0("generating a data frame for plotting for ",i))

		temp_df <- data.frame(age_corr_new=x, y_val=smth_gam$fit + interp_gam[[1]], ci_up=smth_gam$fit + interp_gam[[1]] + orig.ci, ci_dn=smth_gam$fit + interp_gam[[1]] - orig.ci, color_tab=rep(paste0("#",col_tab[j,8]), each=nt), DTI=as.factor(rep(uniDTI[i], each=nt)), region=as.factor(rep(uniReg[j], each=nt)))
		df_all <- rbind(df_all,temp_df)
	}

	df_comb=rbind(df_comb, df_all)
	
    if (ifRelative == 1) 
    {y_label=paste0("Relative ", y_label, sep="")} else if (ifRelative ==2)
    {
		y_label=paste0("Standardized ", y_label, sep="")
	}
    
	nam<-paste("p",i, sep="")
	coln<- y_label
	#attach(df_all)
	assign(
	nam, ggplot() +
	#geom_ribbon(data=df_all, aes_(x=age_corr_new, ymin=ci_dn, ymax=ci_up, group=region), fill=color_tab, alpha=0.3) +
	geom_line(data=df_all, aes_(x=df_all$age_corr_new, y=df_all$y_val, group=df_all$region,  color=df_all$color_tab), size=1, show.legend=FALSE) +
	scale_color_manual(values=c(paste0("#",col_tab[,8]))) +
	labs(x=x_label, y=y_label) +
	xlim(-10, 200) +
	theme(axis.text=element_text(size=16), axis.title = element_text(size = 18)) +
	ylim(ylim_val[[1]], ylim_val[[2]]) 
	)
	#plist<-mget(paste0("p", 1))
	#plist
	df_tmp=data.frame(t(zero_t))
	df_zero<-rbind(df_zero, df_tmp)
	print(paste0("y_label: ", y_label, sep=""))
}


df_dr1=data.frame()
#if (ifRelative == 0)
#{
#=================================
		# plot 1st derivatives (%)
		for (i in 1:length(uniDTI))
		{
			df_all=data.frame()
			rm(zero_t)
			zero_t=array(data=NA, dim=length(uniReg))
			for (j in 1:length(uniReg)) 
			{
				print(paste0(uniDTI[i],' : ',uniReg[j]));

				tmp0=predict(sm.spline(smth_gam_all[[i,j]]$terms, smth_gam_all[[i,j]]$fit+interp_gam_all[[i,j]][[1]], norder=norder_val, df=df_val), x_ext)
				tmp1=predict(sm.spline(smth_gam_all[[i,j]]$terms, smth_gam_all[[i,j]]$fit+interp_gam_all[[i,j]][[1]], norder=norder_val, df=df_val), x_ext,1)
				
				if (cr_per == 0)
				{ y_int_dr1=tmp1
				} else {
				y_int_dr1=100*tmp1/tmp0
				}

				rm(op)
				op=smth.b.all_sp_dr1[i,j,,]    
				tmp <- apply(op, 1, quantile, na.rm=TRUE, probs = c(0.025, 0.975))
				dr1_per.ci=(tmp[2,] - tmp[1,])/2

				x_label=paste("Age (days)")
				xlim_val=c(-10, 200)
				
				if (ifRelative ==0)
				{
					y_label=paste0(toupper(uniDTI[i]))
					if ( i == 1) {
						ylim_val=c(-0.1, 0.7)
					}  else {
						ylim_val=c(-0.35, 0.075)
					}
				
				} else if (ifRelative == 1) {
					y_label=paste0(toupper(uniDTI[i])," C.R. (%)")
					if (as.character(uniDTI[i]) == "tr"){
					y_label="MD C.R (%)"
					}
						ylim_val=c(-0.35, 0.075)
						
					} else {
					y_label=paste0("Standardized ", toupper(uniDTI[i])," C.R. (%)")
					if (as.character(uniDTI[i]) == "tr"){
					y_label="Standardized MD C.R (%)"
					}
						ylim_val=c(-2, 0.2)
				}

				# find the point where curves intersect with zero
				zero_dr1_x=array(data=NA, dim=length(x))
				zero_dr1_y=zero_dr1_x
				ciup_dr1=y_int_dr1 + dr1_per.ci
				cidn_dr1=y_int_dr1 - dr1_per.ci
				for ( t in 1:length(x))
				{
					if ((ciup_dr1[t]>0) & (cidn_dr1[t]<0))
					{
						zero_dr1_x[j]=x[t]
						zero_t[j]=x[t]
						zero_dr1_y[j]=y_int_dr1[t]
						break
					}	
				}
				#generating a data frame for plotting

				temp_df <- data.frame(age_corr_new=x, y_val=y_int_dr1, ci_up_dr1=ciup_dr1, ci_dn_dr1=cidn_dr1, zero_x=zero_dr1_x, zero_y=zero_dr1_y, color_tab=rep(paste0("#",col_tab[j,8]), each=nt), region=as.factor(rep(uniReg[j], each=nt)), DTI=as.factor(rep(uniDTI[i], each=nt)))
				df_all <- rbind(df_all,temp_df)
				
				#write.csv(zero_dr1_x, filename_opt_csv, append<-TRUE)
			}
			df_dr1=rbind(df_dr1, df_all)
			
			nam<-paste("p",i+4, sep="")
			coln<- y_label
			#attach(df_all)
			assign(
			nam, ggplot() +
			#geom_ribbon(data=df_all, aes_(x=df_all$age_corr_new, ymin=df_all$ci_dn_dr1, ymax=df_all$ci_up_dr1, fill=df_all$region), alpha=.3) +
			geom_line(data=df_all, aes_(x=df_all$age_corr_new, group=df_all$region, y=df_all$y_val, color=df_all$color_tab),  size=1, show.legend=FALSE) +
			scale_color_manual(values=c(paste0("#",col_tab[,8]))) +
			geom_point(data=df_all, aes_(x=df_all$zero_x, y=df_all$zero_y,  group=df_all$region), color=df_all$color_tab, na.rm=TRUE, size=5, show.legend=FALSE) +
			labs(x=x_label, y=y_label) +
			xlim(-10, 200) +
			ylim(ylim_val[[1]], ylim_val[[2]]) +
			theme(axis.text=element_text(size=16), axis.title = element_text(size = 18)) +
			geom_hline(yintercept=0, color='black', size=1)
			)
			#plist<-mget(paste0("p", i+4))
			#plist
			
			df_tmp=data.frame(t(zero_t))
			df_zero<-rbind(df_zero, df_tmp)
		}

        row.names(df_zero)=c("FA", "MD", "AD", "RD","dr1% FA", "dr1% MD", "dr1% AD", "dr1% RD")
		colnames(df_zero)=uniReg
		write.csv(df_zero, filename_opt_csv, append<-TRUE, row.names=TRUE, col.names=TRUE)
		
		pdf(file=paste(path,"/",base_name,"_11tracts_gamm_fitMethod",fit_method,"_norder",norder_val,"_gamma",gamma_val,"_m",m_val,"_sp",sp_val,"_k",k_val,"_bs",bs_val,"_df",df_val,"_Relative",ifRelative,"_cr_per",cr_per,"_all.pdf", sep=""), width=20, height=8)
		plist<-mget(paste0("p", 1:8))
		do.call(grid.arrange, c(plist,nrow=2))
		dev.off()
#} else {
	
	#    pdf(file=paste(path,"/",base_name,"_11tracts_gamm_fitMethod",fit_method,"_norder",norder_val,"_gamma",gamma_val,"_m",m_val,"_sp",sp_val,"_k",k_val,"_bs",bs_val,"_df",df_val,"_Relative", ifRelative,"_orig_crPer_only_all.pdf", sep=""), width=20, height=4)
	#	plist<-mget(paste0("p", 1:4))
	#	do.call(grid.arrange, c(plist,nrow=1))
	#	dev.off()
	
#}


print("identify the relationships between absolute original DTI values and its dynamics (original values /per day...")
age_corr_new_uni=unique(df_comb$age_corr_new)
stats_test=data.frame()
for ( i in 1:length(uniDTI))
{
	rm(df_t)
	df_t=data.frame()
	for (t in 1: length(x_ext))
	{
		rm(tmp_orig,tmp_dr1)
		tmp_orig=df_comb[df_comb$DTI==uniDTI[i] & df_comb$age_corr_new==age_corr_new_uni[t],]
		tmp_dr1=df_dr1[df_dr1$DTI==uniDTI[i] & df_dr1$age_corr_new==age_corr_new_uni[t],]
        rm(s_tmp, time, pearson_r, pearson_cidn, pearson_ciup, spearman_r, DTI)
        time = age_corr_new_uni[t]
        pearson_r=cor.test(tmp_orig$y_val, tmp_dr1$y_val, method=c("pearson"))$estimate
        pearson_cidn=cor.test(tmp_orig$y_val, tmp_dr1$y_val, method=c("pearson"))$conf.int[1]
        pearson_ciup=cor.test(tmp_orig$y_val, tmp_dr1$y_val, method=c("pearson"))$conf.int[2]
        spearman_r=cor.test(tmp_orig$y_val, tmp_dr1$y_val, method=c("spearman"))$estimate
        pearson_p=cor.test(tmp_orig$y_val, tmp_dr1$y_val, method=c("pearson"))$p.value
        DTI=uniDTI[i]
        #DTI=uniReg[j]
        s_tmp<-data.frame(time=time, pearson_p=pearson_p, pearson_r=pearson_r, pearson_cidn=pearson_cidn, pearson_ciup=pearson_ciup, spearman_r=spearman_r, DTI=DTI)
        df_t=rbind(df_t, s_tmp)
	}
	stats_test=rbind(stats_test,df_t)
	nam<-paste("p",i, sep="")
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


pdf(file=paste(path,"/",base_name,"_corr_bet_orig_and_speed_orig_Relative",ifRelative,"_cr_per", cr_per,"_all.pdf", sep=""), width=20, height=4)
plist<-mget(paste0("p", 1:4))
do.call(grid.arrange, c(plist,nrow=1))
dev.off()

