library(tcltk2)

win1 <- tktoplevel()
win1$env$rb1 <- tk2radiobutton(win1)
win1$env$rb2 <- tk2radiobutton(win1)
win1$env$rb3 <- tk2radiobutton(win1)

tkbind(win1$env$entName, "<Return>", function() tkdestroy(win1) )
#===============================================
# Define Modflow Cell by Cell Budget Term file
#===============================================
MFmodel = c('ECFTX', 'NPALM', 'LWCSIM')
buttons = c('win1$env$rb1','win1$env$rb2','win1$env$rb3')

rbValue <- tclVar("NPALM")
#for( i in (1:length(MFmodel))){
#  tkconfigure(buttons[i], variable = rbValue, value = MFmodel[i])
#}
tkconfigure(win1$env$rb1, variable = rbValue, value = "NPALM")
tkconfigure(win1$env$rb2, variable = rbValue, value = "ECFTX")
tkconfigure(win1$env$rb3, variable = rbValue, value = "LWCSIM")
  tkgrid(tk2label(win1, text = "Select Model?"),  columnspan = 2, padx = 10, pady = c(15, 5))
  tkgrid(tk2label(win1, text = "NPALM"), win1$env$rb1,  padx = 10, pady = c(0, 5))
  tkgrid(tk2label(win1, text = "ECFTX"), win1$env$rb2,  padx = 10, pady = c(0, 15))
  tkgrid(tk2label(win1, text = "LWCSIM"), win1$env$rb3,  padx = 10, pady = c(0, 15))

win1$env$butOK <- tk2button(win1, text = "OK", width = -6, command = function()  tkdestroy(win1) )

tkgrid(win1$env$butOK, columnspan = 2, padx = 10, pady = c(5, 15))

#tkfocus(win1)
tkraise(win1)
tkwait.window(win1)
MFmodel <-(tclvalue(rbValue))
cbbFile<-switch(MFmodel, 
       "NPALM" = choose.files(default = "\\\\whqhpc01p\\hpcc_shared\\krodberg\\NPB\\*.*"),
       "ECFTX" = choose.files(default = "\\\\whqhpc01p\\hpcc_shared\\dbandara\\CFWI\\ECFTX\\Model\\Transient\\*.*"),
       "LWCSIM" = choose.files(default = "\\\\ad.sfwmd.gov\\dfsroot\\data\\wsd\\MOD\\LWCSIM\\*.*")
       )

