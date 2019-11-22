arquivo = "dados.in"
with open(arquivo, 'r') as arq: #Abre e ja fecha o arquivo r=read
	entrada = arq.readlines() 	#Le todas linhas do arquivo
	maxit = list(entrada[6].split()) #Maxit
	dt = list(entrada[8].split())
	cgd = list(entrada[12].split())
	numArquivos = int(float(maxit[2])/float(cgd[2])) #Converte pra inteiro o numero de execucoes
	
	imax = list(entrada[2].split())
	jmax = list(entrada[4].split())

	xmin =list(entrada[38].split())
	mmin =float(xmin[1][0:len(xmin[1])-2])
	xmax =list(entrada[39].split())
	mmax =float(xmax[1][0:len(xmax[1])-2])

	ymin =list(entrada[41].split())
	nmin =float(ymin[1][0:len(ymin[1])-2])
	ymax =list(entrada[42].split())
	nmax =float(ymax[1][0:len(ymax[1])-2])

	dx = float(mmax-mmin)/(float(imax[2])-1)
	dy = float(nmax-nmin)/(float(jmax[2])-1)
	
	NRBCPoints = list(entrada[73].split())
	NRBC = float(NRBCPoints[2])
	
	PMLx = NRBC*dx
	PMLy = NRBC*dy
	
	sizeX = 1200
	sizeY = 400 #round((nmax-nmin)/(mmax-mmin),1)*sizeX*1.3*2
	yrA = 0 #Ajuste no Yrange

arquivo = "minmax.dat"
with open(arquivo, 'r') as arq2: #Abre e ja fecha o arquivo r=read
	entrada2 = arq2.readlines() 	#Le todas linhas do arquivo
	ListMin = list(entrada2[0].split())
	ListMax = list(entrada2[1].split())
	RhoMin = float(ListMin[0])
	UMin = float(ListMin[1])
	VMin = float(ListMin[2])
	PMin = float(ListMin[3])
	VorMin = float(ListMin[4])
	RhoMax = float(ListMax[0])
	UMax = float(ListMax[1])
	VMax = float(ListMax[2])
	PMax = float(ListMax[3])
	VorMax = float(ListMax[4])
	TMin = float(ListMin[5])
	TMax = float(ListMax[5])
	SOMin = float(ListMin[6])
	SOMax = float(ListMax[6])

min = VorMin #-2
max = VorMax
step = (max-min)/15

nomeArquivo = "plotvorticity.gnu"
arq = open(nomeArquivo, 'w')
arq.write('reset\n')
arq.write('set term jpeg enhanced large font times 15 lw 2 size ' +str(sizeX)+','+str(sizeY)+'\n')
arq.write('#set key samplen 1\n')
arq.write('#set key spacing 1\n')
arq.write('#set cbtics 0.5\n')
#arq.write('set contour base\n')
arq.write('unset sur\n')
#arq.write('set cntrparam levels  incr '+ str(min+step*4) + ',' + str(step) + ',' + str(max-step*4) +' \n')
##arq.write('set cntrparam levels  incr -1.5,' + str((-1.5-0.5)/10) + ',0.5 \n')
arq.write('set view map \n')
arq.write('set iso 100000 \n')
arq.write('set samp 100000 \n')
arq.write("set xrange[0:" + str(mmax-PMLx) + "]\n")
arq.write('set xtics 1 \n')
arq.write('set ytics 0.5 \n')
#arq.write("set yrange["+str(round((nmin+PMLy),2)+yrA)+":"+str(round((nmax-PMLy),2)-yrA)+"] \n")
arq.write("set yrange[-1:1] \n")
#arq.write("set yrange["+str(round((nmin+PMLy),2)+yrA)+":"+str(round((nmax-PMLy),2)-yrA)+"] \n")
arq.write('set zero 1e-12 \n')
arq.write("set size ratio " + str(round((nmax-nmin)/(mmax-mmin),2)) +"  \n")
arq.write('unset key \n')
#arq.write('set palette rgb 22,13,-31 \n')
arq.write('set palette rgb 21,22,23 \n')
arq.write('set pm3d at b \n')
#arq.write('set cbrange['+str(min)+':'+str(max)+'] \n')
arq.write('set cbrange[-1.9:0] \n')

for i in range(0, numArquivos):
	arq.write("set title 'Vorticity contours, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + "' \n")
	arq.write("set output 'vorticity, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + ".jpeg'\n")
	arq.write("splot 'contorno_" + str(i) + ".dat'  u 1:2:7 \n")
arq.close()

min = PMin #0.75
max = PMax #0.65
step = (max-min)/10


nomeArquivo = "plotpressure.gnu"
arq = open(nomeArquivo, 'w')
arq.write('reset\n')
arq.write('set term jpeg enhanced large font times 20 lw 2 size ' +str(sizeX)+','+str(sizeY)+'\n')
arq.write('#set key samplen 1\n')
arq.write('#set key spacing 1\n')
arq.write('#set cbtics 0.5\n')
#arq.write('set contour base\n')
arq.write('unset sur\n')
#arq.write('set cntrparam levels  incr '+ str(min) + ',' + str((max-min)/10) + ',' + str(max) +' \n')
arq.write('set view map \n')
arq.write('set iso 100000 \n')
arq.write('set samp 1000000 \n')
arq.write("set xrange[0:" + str(mmax-PMLx) + "]\n")
arq.write('set xtics 1 \n')
arq.write('set ytics 0.5 \n')
#arq.write("set yrange["+str(round((nmin+PMLy),2)+yrA)+":"+str(round((nmax-PMLy),2)-yrA)+"] \n")
arq.write("set yrange[-1:1] \n")
arq.write('set zero 1e-12 \n')
arq.write("set size ratio " + str(round((nmax-nmin)/(mmax-mmin),2)) +"  \n")
arq.write('unset key \n')
arq.write('set palette rgb 33,13,10 \n')
arq.write('set pm3d at b \n')
#arq.write('set cbrange[-1.9:0] \n')
arq.write('set cbrange['+str(min)+':'+str(max)+'] \n')

for i in range(0, numArquivos):
	arq.write("set title 'Pressure contours, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + "' \n")
	arq.write("set output 'pressure, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + ".jpeg'\n")
	arq.write("splot 'contorno_" + str(i) + ".dat'  u 1:2:6 \n")
arq.close()

min = VMin #-0.05
max = VMax #0.05


nomeArquivo = "plotyvel.gnu"
arq = open(nomeArquivo, 'w')
arq.write('reset\n')
arq.write('set term jpeg enhanced large font times 20 lw 2 size ' +str(sizeX)+','+str(sizeY)+'\n')
arq.write('#set cbtics 0.5\n')
arq.write('#set key samplen 1\n')
arq.write('#set key spacing 1\n')
#arq.write('set contour base\n')
arq.write('unset sur\n')
#arq.write('set cntrparam levels  incr '+ str(min) + ',' + str((max-min)/10) + ',' + str(max) +' \n')
##arq.write('set cntrparam levels  incr -0.04,  ,0.04 \n')
arq.write('set view map \n')
arq.write('set iso 100000 \n')
arq.write('set samp 1000000 \n')
arq.write("set xrange[0:" + str(mmax-PMLx) + "]\n")
arq.write('set xtics 1 \n')
arq.write('set ytics 0.5 \n')
arq.write("set yrange[-1:1] \n")
#arq.write("set yrange["+str(round((nmin+PMLy),2)+yrA)+":"+str(round((nmax-PMLy),2)-yrA)+"] \n")
arq.write('set zero 1e-12 \n')
arq.write("set size ratio " + str(round((nmax-nmin)/(mmax-mmin),2)) +"  \n")
arq.write('unset key \n')
arq.write('set palette rgb 33,13,10 \n')
arq.write('set pm3d at b \n')
arq.write('set cbrange['+str(min)+':'+str(max)+'] \n')
#arq.write('set cbrange[-1.9:0] \n')

for i in range(0, numArquivos):
	arq.write("set title 'Transversal velocity contours, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + "' \n")
	arq.write("set output 'Transversal velocity, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + ".jpeg'\n")
	arq.write("splot 'contorno_" + str(i) + ".dat'  u 1:2:5 \n")
arq.close()

min = UMin
max = UMax


nomeArquivo = "plotxvel.gnu"
arq = open(nomeArquivo, 'w')
arq.write('reset\n')
arq.write('set term jpeg enhanced large font times 20 lw 2 size ' +str(sizeX)+','+str(sizeY)+'\n')
arq.write('#set key samplen 1\n')
arq.write('#set key spacing 1\n')
#arq.write('set contour base\n')
arq.write('unset sur\n')
#arq.write('set cntrparam levels  incr '+ str(min) + ',' + str((max-min)/10) + ',' + str(max) +' \n')
##arq.write('set cntrparam levels  incr -0.5, 0.05 ,1.5 \n')
arq.write('set view map \n')
arq.write('set iso 100000 \n')
arq.write('set samp 1000000 \n')
arq.write("set xrange[0:" + str(mmax-PMLx) + "]\n")
arq.write('set xtics 1 \n')
arq.write('set ytics 0.5 \n')
arq.write("set yrange[-1:1] \n")
#arq.write("set yrange["+str(round((nmin+PMLy),2)+yrA)+":"+str(round((nmax-PMLy),2)-yrA)+"] \n")
arq.write('set zero 1e-12 \n')
arq.write("set size ratio " + str(round((nmax-nmin)/(mmax-mmin),2)) +"  \n")
arq.write('unset key \n')
arq.write('set palette rgb 33,13,10 \n')
arq.write('set pm3d at b \n')
arq.write('set cbrange['+str(min)+':'+str(max)+'] \n')
#arq.write('set cbrange[-0.5:1.5] \n')

for i in range(0, numArquivos):
	arq.write("set title 'Longitudinal velocity contours, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + "' \n")
	arq.write("set output 'Longitudinal velocity, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + ".jpeg'\n")
	arq.write("splot 'contorno_" + str(i) + ".dat'  u 1:2:4 \n")
arq.close()


min = RhoMin
max = RhoMax
step = (max-min)/15


nomeArquivo = "plotrho.gnu"
arq = open(nomeArquivo, 'w')
arq.write('reset\n')
arq.write('set term jpeg enhanced large font times 20 lw 2 size ' +str(sizeX)+','+str(sizeY)+'\n')
arq.write('#set key samplen 1\n')
arq.write('#set key spacing 1\n')
#arq.write('set contour base\n')
arq.write('unset sur\n')
#arq.write('set cntrparam levels  incr '+ str(min+step) + ',' + str((max-min)/10) + ',' + str(max-step) +' \n')
##arq.write('set cntrparam levels  incr 1,' + str((1.3-1)/10) +',1.3 \n')
arq.write('set view map \n')
arq.write('set iso 100000 \n')
arq.write('set samp 1000000 \n')
arq.write("set xrange[0:" + str(mmax-PMLx) + "]\n")
arq.write('set xtics 1 \n')
arq.write('set ytics 0.5 \n')
arq.write("set yrange[-1:1] \n")
#arq.write("set yrange["+str(round((nmin+PMLy),2)+yrA)+":"+str(round((nmax-PMLy),2)-yrA)+"] \n")
arq.write('set zero 1e-12 \n')
arq.write("set size ratio " + str(round((nmax-nmin)/(mmax-mmin),2)) +"  \n")
arq.write('unset key \n')
arq.write('set palette rgb 22,13,-31 \n')
arq.write('set pm3d at b \n')
arq.write('set cbrange['+str(min)+':'+str(max)+'] \n')
#arq.write('set cbrange[1:1.3] \n')

for i in range(0, numArquivos):
	arq.write("set title 'Specific Mass contours, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + "' \n")
	arq.write("set output 'Specific Mass, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + ".jpeg'\n")
	arq.write("splot 'contorno_" + str(i) + ".dat'  u 1:2:3 \n")
arq.close()

min = TMin
max = TMax
step = (max-min)/12

nomeArquivo = "plottemp.gnu"
arq = open(nomeArquivo, 'w')
arq.write('reset\n')
arq.write('set term jpeg enhanced large font times 20 lw 2 size ' +str(sizeX)+','+str(sizeY)+'\n')
arq.write('#set key samplen 1\n')
arq.write('#set key spacing 1\n')
#arq.write('set contour base\n')
arq.write('unset sur\n')
#arq.write('set cntrparam levels  incr '+ str(min+step*4) + ',' + str(step) + ',' + str(max-step*4) +' \n')
##arq.write('set cntrparam levels  incr 0.5, ' + str((0.7-0.5)/10) + ' ,0.7 \n')
arq.write('set view map \n')
arq.write('set iso 100000 \n')
arq.write('set samp 1000000 \n')
arq.write("set xrange[0:" + str(mmax-PMLx) + "]\n")
arq.write('set xtics 1 \n')
arq.write("set yrange[-1:1] \n")
arq.write('set ytics 0.5 \n')
arq.write("set yrange[-1:1] \n")
#arq.write("set yrange["+str(round((nmin+PMLy),2)+yrA)+":"+str(round((nmax-PMLy),2)-yrA)+"] \n")
arq.write('set zero 1e-12 \n')
arq.write("set size ratio " + str(round((nmax-nmin)/(mmax-mmin),2)) +"  \n")
arq.write('unset key \n')
arq.write('set palette rgb 22,13,-31 \n')
arq.write('set pm3d at b \n')
arq.write('set cbrange['+str(min)+':'+str(max)+'] \n')

for i in range(0, numArquivos):
	arq.write("set title 'Temperature contours, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + "' \n")
	arq.write("set output 'Temperature, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + ".jpeg'\n")
	arq.write("splot 'contorno_" + str(i) + ".dat'  u ($1):($2):(1.4*($6)/($3))\n")
arq.close()

min = SOMin
max = SOMax
step = (max-min)/12

nomeArquivo = "plotO2.gnu"
arq = open(nomeArquivo, 'w')
arq.write('reset\n')
arq.write('set term jpeg enhanced large font times 20 lw 2 size ' +str(sizeX)+','+str(sizeY)+'\n')
arq.write('#set key samplen 1\n')
arq.write('#set key spacing 1\n')
#arq.write('set contour base\n')
arq.write('unset sur\n')
#arq.write('set cntrparam levels  incr '+ str(min+step*4) + ',' + str(step) + ',' + str(max-step*4) +' \n')
##arq.write('set cntrparam levels  incr 0.5, ' + str((0.7-0.5)/10) + ' ,0.7 \n')
arq.write('set view map \n')
arq.write('set iso 100000 \n')
arq.write('set samp 1000000 \n')
arq.write("set xrange[0:" + str(mmax-PMLx) + "]\n")
arq.write('set xtics 1 \n')
arq.write("set yrange[-1:1] \n")
arq.write('set ytics 0.5 \n')
arq.write("set yrange[-1:1] \n")
#arq.write("set yrange["+str(round((nmin+PMLy),2)+yrA)+":"+str(round((nmax-PMLy),2)-yrA)+"] \n")
arq.write('set zero 1e-12 \n')
arq.write("set size ratio " + str(round((nmax-nmin)/(mmax-mmin),2)) +"  \n")
arq.write('unset key \n')
arq.write('set palette rgb 22,13,-31 \n')
arq.write('set pm3d at b \n')
arq.write('set cbrange['+str(min)+':'+str(max)+'] \n')

for i in range(0, numArquivos):
	arq.write("set title 'Oxigen contours, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + "' \n")
	arq.write("set output 'Oxigen, t=" + str(int(float(cgd[2])*float(dt[2][1:len(dt[2])-2]))*(i+1)) + ".jpeg'\n")
	arq.write("splot 'contorno_" + str(i) + ".dat'  u ($1):($2):($8)\n")
arq.close()

