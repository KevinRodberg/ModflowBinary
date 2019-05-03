win1 <- tktoplevel()
MFmodel = c('ECFTX', 'NPALM', 'LWCSIM')

numIDs = length(MFmodel);
bs = vector();
tkgrid(tk2label(win1, text = "Select Model?"),  columnspan = 2, padx = 10, pady = c(15, 5))
for (num in 1:numIDs) {
  b <- tk2radiobutton(win1);
  tkconfigure(b, variable = rbValue, value = MFmodel[num])
  tkgrid(tk2label(win1, text = MFmodel[num]), b,  padx = 10, pady = c(0, 5))
  bs= append(bs,b)
}
butOK <- tk2button(win1, text = "OK", width = -6, command = function()  tkdestroy(win1) )
tkgrid(butOK, columnspan = 2, padx = 10, pady = c(5, 15))

tkraise(win1)
tkwait.window(win1)
model  <-(tclvalue(rbValue))
cbbFile<-switch(model, 
                "NPALM" = choose.files(default = "\\\\whqhpc01p\\hpcc_shared\\krodberg\\NPB\\*.*"),
                "ECFTX" = choose.files(default = "\\\\whqhpc01p\\hpcc_shared\\dbandara\\CFWI\\ECFTX\\Model\\Transient\\*.*"),
                "LWCSIM" = choose.files(default = "\\\\ad.sfwmd.gov\\dfsroot\\data\\wsd\\MOD\\LWCSIM\\*.*")
)
