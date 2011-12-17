#!/usr/bin/python
import numpy
import pylab
import scipy.special
import math
import matplotlib.tri as tri

coeff=32
omega=1.0/(0.5+0.5/coeff)
radius=39
diffusion=1.0/3.0*(1.0/omega-0.5)
num_terms=20

def read_file(file_name):
    concentration=numpy.loadtxt(file_name)
    pylab.imshow(concentration)
    

def read_files(file_dir):
    styles=["kv","k^","ks","ko","kD"]    
    non_time=[]    
    legs=[]    
    for num,counter in enumerate(range(200*coeff,2000*coeff,400*coeff)):
        file_name=file_dir+"density"+(7-len(str(counter)))*"0"+str(counter)+".dat"
        concentration=numpy.loadtxt(file_name)
        dims=concentration.shape
        #pylab.figure()
        #pylab.imshow(concentration)
        x=numpy.arange(0,dims[1])        
        x_rad=numpy.linspace(0.0,1.0,radius+1)

        
        bessel=numpy.zeros_like(x)
        bessel_zeros=scipy.special.jn_zeros(0,num_terms)
        non_time.append(diffusion*counter/(radius*radius))
        for term in range(0,num_terms):
            bessel=bessel+2/(bessel_zeros[term]*scipy.special.jn(1,bessel_zeros[term]))\
                *numpy.exp(-bessel_zeros[term]*bessel_zeros[term]*diffusion*counter/(radius*radius))\
                *scipy.special.jn(0,bessel_zeros[term]*x/radius)
        bessel=1-bessel
        #bessel[numpy.where(bessel)>1]=0
        #bessel=numpy.roll(bessel,dims[1]/2)
        print "Shape",dims        
        pylab.figure(1)
        pylab.plot(x_rad,concentration[(dims[0]-1)/2,(dims[1]-1)/2:(dims[1]-1)/2+radius+1],styles[num])        
        #pylab.plot(concentration[dims[0]/2-1,dims[1]/2:dims[1]/2+radius])
        #pylab.figure(2)
        pylab.plot(x_rad,bessel[numpy.where(x<=radius)],styles[num]+"--",markerfacecolor="None")
    legs=[[r'''$\tau='''+str(time)[0:6]+'''$''',r'''$\tau='''+str(time)[0:6]+'''$'''] for time in non_time]    
    legs=numpy.ravel(legs)   
    pylab.title(r'''$D='''+str(diffusion)[0:6]+'''$''',fontsize=30)
    pylab.legend(legs,fancybox=True,loc=2)
    pylab.xlabel(r'''$x$''',fontsize=20)
    pylab.ylabel(r'''$C$''',fontsize=20)
    #pylab.legend(legs,fancybox=True,labelspacing=0.1)
    pylab.savefig("cylinder"+str(diffusion)[2:6]+".eps",dpi=300)        
        
def show_bessel():
    x=numpy.arange(0.0,10.0,0.1)
    bessel=scipy.special.jn(0,x)
    pylab.plot(x,bessel)

def read_film(file_dir):
    concentration=numpy.loadtxt(file_dir+"conc_initial.dat")
    ux=numpy.loadtxt(file_dir+"ux_initial.dat")
    uy=numpy.loadtxt(file_dir+"uy_initial.dat")
    print concentration.shape
    pylab.figure()    
    pylab.imshow(concentration)
    pylab.figure()
    pylab.plot(concentration[:,200])
    pylab.figure()
    pylab.plot(concentration[:,0])
    pylab.figure()
    pylab.imshow(uy)
    pylab.colorbar()

def film_analytical():
    num_terms=30
    c0=0.0
    cs=1.0
    u_bubble=0.05
    y,x=0.01*numpy.mgrid[1:101,1:501]
    
    c=numpy.zeros_like(x)
    diffusion=1.0/3.0*(1-0.5)
    
    sign=-1
    for i in range(0,num_terms):
        sign=-sign        
        c=c+sign*(scipy.special.erfc((y+2*i)/numpy.sqrt(4*x*diffusion/u_bubble))+scipy.special.erfc((2*(i+1)-y)/numpy.sqrt(4*x*diffusion/u_bubble)))
    print "X=",x
    print "Y=",y
    print "C=",c
    c=c0-(c0-cs)*c
    pylab.imshow(c)
    pylab.colorbar()

def compare_film(file_dir):
    concentration=numpy.loadtxt(file_dir+"film0050000.dat")
    uy=numpy.loadtxt(file_dir+"uy_initial.dat")
    
    dims=concentration.shape    
    print dims
    pylab.figure()    
    pylab.imshow(concentration[:,:-4],extent=(0,10,0,1))
    pylab.title("Simulations")
    pylab.colorbar()    
    #pylab.figure()
    #pylab.plot(concentration[:,dims[1]/2])
    #pylab.figure()
    #pylab.plot(concentration[:,0])

    num_terms=30
    c0=0.5
    cs=1.0
    u_bubble=0.05*20   
    y,x=0.01*numpy.mgrid[1:101,1:1001]
    
    c=numpy.zeros_like(x)
    diffusion=1.0/3.0*(1.0/1.4-0.5)
    
    res=y/numpy.sqrt(4*x*diffusion/u_bubble)
    print res.shape
    sign=-1
    for i in range(0,num_terms):
        sign=-sign
        c=c+sign*(scipy.special.erfc((y+2.0*i)/numpy.sqrt(4.0*x*diffusion/u_bubble))+scipy.special.erfc((2.0*(i+1)-y)/numpy.sqrt(4.0*x*diffusion/u_bubble)))
        
    #print "X=",x
    #print "Y=",y
    #print "C=",c
    c=c0-(c0-cs)*c
    pylab.figure()    
    pylab.imshow(c,extent=(0,10,0,1))
    pylab.colorbar()
    pylab.title("Analytics")

    pylab.figure()
    pylab.imshow(uy)
    pylab.colorbar()
    
    pylab.figure()
    pylab.plot(concentration[0,:])
    pylab.figure()
    pylab.plot(concentration[:,0])
    pylab.figure()
    pylab.plot(concentration[:,dims[1]-1])
    
    c_levels=numpy.arange(0.0,1.0,0.1)
    print c_levels
    pylab.figure(figsize=(10,1))
    c1=pylab.contour(concentration[:,:-4],levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c1,fontsize=9, inline=1)    
    #pylab.figure(figsize=(10,1))
    c2=pylab.contour(c,levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c2,fontsize=9, inline=1)    
    
    print numpy.where(concentration>1.8)

def compare_three_films(file_dir):
    conc_antibb =numpy.loadtxt(file_dir+"film0050000.dat")
    conc_inamuro=numpy.loadtxt(file_dir+"film_inamuro0050000.dat")
    conc_outflow=numpy.loadtxt(file_dir+"film_outflow0050000.dat")    
    
    pylab.figure()    
    pylab.imshow(conc_antibb[:,:-20],extent=(0,10,0,1))
    pylab.title("AntiBB simulations")
    pylab.colorbar()    
    pylab.figure()
    pylab.imshow(conc_inamuro[:,:-20],extent=(0,10,0,1))
    pylab.title("Inamuro simulations")
    pylab.colorbar()    
    pylab.figure()
    pylab.imshow(conc_outflow[:,:-20],extent=(0,10,0,1))
    pylab.title("Partial Inamuro simulations")
    pylab.colorbar()

    c_levels=numpy.arange(0.6,1.0,0.1)
    print c_levels
    pylab.figure(99,figsize=(10,1))
    c1=pylab.contour(conc_antibb[:,:-20],levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c1,fontsize=9, inline=1)    
    #pylab.figure(figsize=(10,1))
    c2=pylab.contour(conc_inamuro[:,:-20],levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c2,fontsize=9, inline=1)
    c3=pylab.contour(conc_outflow[:,:-20],levels=c_levels,extent=(0.0,10.0,0.0,1.0))



    num_terms=30
    c0=0.5
    cs=1.0
    u_bubble=0.05*20   
    y,x=0.01*numpy.mgrid[1:101,1:1001]
    
    c=numpy.zeros_like(x)
    diffusion=1.0/3.0*(1.0/1.4-0.5)
    
    res=y/numpy.sqrt(4*x*diffusion/u_bubble)
    print res.shape
    sign=-1
    for i in range(0,num_terms):
        sign=-sign
        c=c+sign*(scipy.special.erfc((y+2.0*i)/numpy.sqrt(4.0*x*diffusion/u_bubble))+scipy.special.erfc((2.0*(i+1)-y)/numpy.sqrt(4.0*x*diffusion/u_bubble)))
        
    c=c0-(c0-cs)*c
    pylab.figure()    
    pylab.imshow(c,extent=(0,10,0,1))
    pylab.colorbar()
    pylab.title("Analytics")

    pylab.figure(99,)
    c_anal=pylab.contour(c[:,:-20],levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c_anal,fontsize=9, inline=1)    


def compare_full_profiles(file_dir):
    conc_antibb =numpy.loadtxt(file_dir+"FullProfile/"+"film_antibb0050000.dat")
    conc_inamuro=numpy.loadtxt(file_dir+"FullProfile/"+"film_outflow0050000.dat")
    
    print "Shape=",conc_antibb.shape
    pylab.figure()    
    pylab.imshow(conc_antibb,extent=(0,20,0,1))
    pylab.title("AntiBB simulations")
    pylab.colorbar()    
    pylab.figure()
    pylab.imshow(conc_inamuro,extent=(0,20,0,1))
    pylab.title("Inamuro simulations")
    pylab.colorbar()    

    c_levels=numpy.arange(0.2,1.0,0.1)
    print c_levels
    pylab.figure(99,figsize=(20,1))
    c1=pylab.contour(conc_antibb,levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c1,fontsize=9, inline=1)    
    #pylab.figure(figsize=(10,1))
    c2=pylab.contour(conc_inamuro,levels=c_levels,extent=(0.0,10.0,0.0,1.0))
    pylab.clabel(c2,fontsize=9, inline=1)



    num_terms=100
    c0=0.0
    cs=1.0
    ny=40
    omega=1.8
    u_bubble=0.05
    diffusion=1.0/3.0*(1.0/omega-0.5)
    pe=u_bubble*ny/diffusion
    print "Pe=",pe
    y,x=0.005*numpy.mgrid[1:201,1:8001]
    #pe=pe/2
    c=numpy.zeros_like(x)
    
    sign=-1
    for i in range(0,num_terms):
        sign=-sign
        c=c+sign*(scipy.special.erfc((y+2.0*i)/numpy.sqrt(4.0*x/pe))+scipy.special.erfc((2.0*(i+1)-y)/numpy.sqrt(4.0*x/pe)))
    c=c0-(c0-cs)*c
    #c=scipy.special.erfc(y/numpy.sqrt(4.0*x/pe))    
    print c    
    pylab.figure()    
    pylab.imshow(c,extent=(0,20,0,0.5))
    pylab.colorbar()
    pylab.title("Analytics")

    pylab.figure(99)
    #pylab.figure(figsize=(10,1))
    c_anal=pylab.contour(c,levels=c_levels,linestyles="dashed",extent=(0.0,10.0,0.0,0.5))
    #pylab.clabel(c_anal,fontsize=9, inline=1)

def film_reconstruction(file_dir):
    #bessel_zeros=scipy.special.jv(-3.0/4.0,num_terms)
    wmfivetwo=numpy.array([0.20996479709362217, 0.018180897590318785, 0.006919303361557072,\
                           0.0037318500468507265, 0.0023673658890243144, \
                           0.0016504522182779553, 0.0012240656323313512, \
                           0.0009483894730813764, 0.0007590968780885121, \
                           0.0006230691899630487, 0.000521773141559365, \
                           0.0004441474396258039, 0.00038324243624026397, \
                           0.00033450473194280715, 0.0002948449039830018, \
                           0.0002621039555673866, 0.0002347343932804438, \
                           0.00021160239472918784, 0.00019186109713872825, \
                           0.00017486713895424112, 0.00016012432055124612, \
                           0.00014724473213250115, 0.0001359214046973261, \
                           0.00012590872750600792, 0.00011700820222555297, \
                           0.0001090579288501799, 0.00010192474302362525, \
                           0.00009549826478267914, 0.00008968634379190777, \
                           0.00008441153749062422, 0.00007960836196536182, \
                           0.00007522112702512196, 0.00007120221729922205, \
                           0.00006751071698186362, 0.0000641113016165315, \
                           0.000060973339051508034, 0.00005807015547318643, \
                           0.00005537843263812609, 0.00005287771007323228, \
                           0.00005054997178480274, 0.000048379301408858354, \
                           0.00004635159310216252, 0.00004445430807096168, \
                           0.00004267626865559752, 0.00004100748346843951, \
                           0.00003943899832564699, 0.000037962768698096496, \
                           0.000036571550189571094, 0.00003525880417636117, \
                           0.00003401861624793521, 0.00003284562549298673, \
                           0.00003173496300815666, 0.000030682198273720266, \
                           0.00002968329226266342, 0.0000287345563299556, \
                           0.000027832616078876215, 0.0000269743795242561, \
                           0.00002615700897675555, 0.00002537789615694936, \
                           0.000024634640120406708, 0.00002392502763531586, \
                           0.000023247015704908107, 0.000022598715969558762, \
                           0.00002197838076039508, 0.000021384390606644113, \
                           0.000020815243025254824, 0.00002026954244416346, \
                           0.000019745991129333714, 0.00001924338100264011, \
                           0.000018760586251353425, 0.000018296556642761726, \
                           0.000017850311467667332, 0.000017420934045730063, \
                           0.000017007566733802623, 0.000016609406384810287, \
                           0.000016225700211394876, 0.0000158557420132339, \
                           0.000015498868731945504, 0.000015154457301158608, \
                           0.000014821921763302157, 0.00001450071062731367, \
                           0.000014190304444456835, 0.000013890213582018029, \
                           0.000013599976176311622, 0.000013319156248667392, \
                           0.000013047341969944055, 0.000012784144059735866, \
                           0.00001252919430900766, 0.000012282144214937344, \
                           0.000012042663718404304, 0.000011810440035547778, \
                           0.00001158517657507229, 0.000011366591934418066, \
                           0.00001115441896833911, 0.000010948403923597011, \
                           0.000010748305634853729, 0.000010553894776504868, \
                           0.000010364953166337604, 0.000010181273116649793, \
                           0.000010002656829307337, 9.828915831386193e-6, \
                           9.659870448208707e-6, 9.49534931097945e-6, \
                           9.335188896470898e-6, 9.179233096398142e-6, \
                           9.027332814266679e-6, 8.879345587644758e-6, \
                           8.73513523417869e-6, 8.59457151961857e-6, 8.457529846050854e-6, \
                           8.32389095943517e-6, 8.1935406745277e-6, 8.066369616458602e-6, \
                           7.942272977545397e-6, 7.821150288441995e-6, \
                           7.702905202654277e-6, 7.5874452934925025e-6, \
                           7.47468186263942e-6, 7.364529759619492e-6, \
                           7.256907211476081e-6, 7.151735661952972e-6, \
                           7.048939619510668e-6, 6.94844651384123e-6, \
                           6.850186560093698e-6, 6.754092630545741e-6, \
                           6.660100132994738e-6, 6.568146895868383e-6, \
                           6.478173059186307e-6, 6.390120971368812e-6, \
                           6.3039350913559065e-6, 6.219561895780246e-6, \
                           6.136949790954026e-6, 6.056049029226355e-6, \
                           5.976811629690077e-6, 5.89919130273287e-6, \
                           5.823143378440544e-6, 5.748624738477849e-6, \
                           5.675593751315735e-6, 5.604010210652239e-6, \
                           5.533835276712681e-6, 5.465031420482145e-6, \
                           5.397562370506302e-6, 5.331393062241362e-6, \
                           5.266489589791824e-6, 5.202819159859141e-6, \
                           5.140350047878658e-6, 5.079051556127257e-6, \
                           5.01889397376408e-6, 4.959848538681422e-6, \
                           4.901887401087894e-6, 4.844983588703604e-6, \
                           4.789110973507299e-6, 4.734244239954024e-6, \
                           4.680358854602334e-6, 4.627431037040406e-6, \
                           4.575437732101662e-6, 4.524356583247067e-6, \
                           4.4741659071123195e-6, 4.424844669102156e-6, \
                           4.376372460064756e-6])
    wmcubic=numpy.array([0.13797406599320738, 0.036854581143239175, 0.02133155594370827, \
                         0.015010686087355087, 0.011579602879061325, 0.009425240459319681, \
                         0.00794676552988761, 0.006869235708429393, 0.006049027438296182, \
                         0.005403797343711836, 0.004882949141705082, 0.0044536787382349384, \
                         0.0040937854185415555, 0.0037877079495512375, \
                         0.0035242151877391964, 0.003294997916730558, \
                         0.0030937766542952425, 0.002915717555439377, \
                         0.0027570390413366227, 0.0026147402465559653, \
                         0.0024864094328322945, 0.002370086181472648, 0.002264160539779118, \
                         0.0021672980550348627, 0.0020783832616715044, \
                         0.0019964765310819576, 0.0019207807375925753, \
                         0.001850615230525651, 0.001785395309967969, 0.0017246158947239064, \
                         0.0016678384163521522, 0.0016146802195225812, \
                         0.0015648059267859933, 0.0015179203557653368, \
                         0.0014737626725206412, 0.0014321015367415484, \
                         0.0013927310476398757, 0.0013554673405968029, 0.00132014571598235, \
                         0.0012866182056656857, 0.0012547515008308823, \
                         0.0012244251805145821, 0.0011955301909566858, \
                         0.0011679675352164775, 0.001141647140020983, \
                         0.0011164868724411468, 0.001092411683789792, \
                         0.0010693528619443273, 0.0010472473764049522, \
                         0.0010260373029400856, 0.001005669316760948, \
                         0.0009860942448871525, 0.0009672666697926882, \
                         0.0009491445776069171, 0.0009316890451347559, \
                         0.0009148639607888969, 0.000898635775223389, \
                         0.0008829732780440975, 0.0008678473974712466, \
                         0.0008532310202445551, 0.000839098829416329, \
                         0.0008254271580562615, 0.0008121938569205483, \
                         0.0007993781747559135, 0.0007869606497074265, \
                         0.0007749230107130705, 0.0007632480877511553, \
                         0.0007519197302406965, 0.0007409227323453307, \
                         0.0007302427650298862, 0.0007198663136691806, \
                         0.0007097806209968215, 0.0006999736348672356, \
                         0.0006904339601143648, 0.0006811508144305624, \
                         0.0006721139877036992, 0.0006633138045573109, \
                         0.0006547410897960503, 0.0006463871364960805, \
                         0.0006382436765076748, 0.0006303028531598466, \
                         0.0006225571959772052, 0.0006149995972390405, 0.00060762329023299, \
                         0.0006004218290146798, 0.0005933890696697541, \
                         0.0005865191527910231, 0.0005798064872363093, \
                         0.0005732457349049316, 0.0005668317966193273, \
                         0.0005605597988617798, 0.000554425081500312, 0.000548423186165646, \
                         0.0005425498454906621, 0.0005368009729481291, \
                         0.0005311726534411112, 0.0005256611343459384, \
                         0.0005202628171764819, 0.0005149742497774396, \
                         0.000509792118961795, 0.0005047132435614256, \
                         0.0004997345679090261, 0.0004948531557746551, \
                         0.0004900661844907577, 0.00048537093960920554, \
                         0.00048076480969575636, 0.0004762452815198213, \
                         0.00047180993548151946, 0.0004674564412607595, \
                         0.00046318255379765045, 0.0004589861093689904, \
                         0.00045486502196710586, 0.00045081727983601596, \
                         0.0004468409421967647, 0.00044293413614472334, \
                         0.00043909505370775154, 0.00043532194905709584, \
                         0.0004316131358594908, 0.0004279669847644897, \
                         0.00042438192101768214, 0.00042085642219294103, \
                         0.000417389016036788, 0.00041397827841843323, \
                         0.00041062283137950596, 0.0004073213412778911, \
                         0.0004040725170212972, 0.0004008751083791981, \
                         0.00039772790438793656, 0.0003946297318123248, \
                         0.000391579453691969, 0.000388575967949311, 0.0003856182060619201, \
                         0.000382705131795076, 0.0003798357399915421, \
                         0.00037700905541481736, 0.0003742241316443967, \
                         0.00037148005001910067, 0.0003687759186267475, \
                         0.0003661108713376293, 0.0003634840668795648, \
                         0.0003608946879525167, 0.00035834194038075984, \
                         0.0003558250523007927, 0.00035334327338323315, \
                         0.0003508958740870435, 0.00034848214494454037, \
                         0.0003461013958757237, 0.0003437529555304854, \
                         0.00034143617065746486, 0.0003391504054982293, \
                         0.00033689504120561653, 0.0003346694752851603, \
                         0.0003324731210584933, 0.000330305407147744, \
                         0.00032816577697996523, 0.00032605368831072697, \
                         0.0003239686127659561, 0.0003219100354012626, \
                         0.00031987745427795427, 0.00031787038005491955])
    
    bessel_zeros=numpy.array([1.45499708549819386722507026854180579346, \
                              2.92713300440027766553121136094888110735, \
                              3.8575781013492938200342434118812546441, \
                              4.60177773228129985324241798621793597148, \
                              5.24082406699070314354952418816449605517, \
                              5.80978786426398906487468774457944491812, \
                              6.32769430108396681532057061373337375445, \
                              6.8062480022514212669586769504065468988, \
                              7.25326168965686308822824791463505009483, \
                              7.6742592995607184149250788121278225534, \
                              8.07331795319108946335626947213305371378, \
                              8.45354897437340233163926840148873400306, \
                              8.81739086884080680472896916336389339722, \
                              9.16679702363228751221979864959553140746, \
                              9.50336098897760887069148869968298949199, \
                              9.82840298647447126437751319761265707777, \
                              10.14303137712219911244514916412671544951, \
                              10.44818741811621983414432659838862125807, \
                              10.74467854808127573053205335964669225223, \
                              11.03320360258297538892991564331548075855, \
                              11.31437222998145142057443770168774860219, \
                              11.58872005940677119103217820431909528765, \
                              11.85672070452047449321754322768063874146, \
                              12.11879537438111167446084673341581712359, \
                              12.37532064988645588833414500384166974648, \
                              12.62663483645073127192279608364483685065, \
                              12.8730431991488418101183117224826294002, \
                              13.11482231162738863620077232279337810946, \
                              13.35222369554359348133575969794532279845, \
                              13.58547688707664934803393403887138075079, \
                              13.81479203704230978691485194517723996389, \
                              14.0403621284925403691860805796792028747, \
                              14.2623648784142054972133104428159388272, \
                              14.4809643768493324602515059590286238574, \
                              14.69631250643745383397730190280670992963, \
                              14.90855017729768053412230968849290743382, \
                              15.11780840578937570988470787292523813834, \
                              15.32420926061949995243758582398480085884, \
                              15.52786669570596600871656820859722187773, \
                              15.7288872859366536654040406731452138516, \
                              15.9273708793136173631920062299682027691, \
                              16.12341117681167755300451799359118160072, \
                              16.31709624950989469369530608919286116217, \
                              16.50850900109559236347854935232210542604, \
                              16.69772758263279807557262678658274113287, \
                              16.88482576548234240189867071924417329888, \
                              17.0698732774215122579902242575503858456, \
                              17.25293610630693319483410269200861091738, \
                              17.43407677503112073131534139072114313752, \
                              17.61335459102147480760264415427211297139, \
                              17.79082587310469749719229172428528076748, \
                              17.96654415819694343585843588062836815013, \
                              18.14056038997006420349764940438806368648, \
                              18.31292309137856038366383535925968134061, \
                              18.48367852270329730273712784271426510892, \
                              18.65287082657088177261838964015053009292, \
                              18.82054216123703645850384614058496825444, \
                              18.98673282327435071906554950975388778683, \
                              19.15148136067609592105538219021147335869, \
                              19.31482467727557424266429445652771548091, \
                              19.47679812928237338488715497506428607718, \
                              19.63743561465094432179652946866365109956, \
                              19.79676965592142972660408702565954850112, \
                              19.95483147710622580687699207508067817164, \
                              20.1116510751371510811081913728251719114, \
                              20.26725728633629104179428517626838796877, \
                              20.42167784832770592539536914201674076594, \
                              20.5749394577664739838553947225769077593, \
                              20.72706782422534501022699390743200749406, \
                              20.87808772054703867991961754439254908311, \
                              21.02802302994145620306095998066847616749, \
                              21.17689679008136355778581681928277317273, \
                              21.32473123442708740654118303228698174865, \
                              21.47154783099012537957908642158269803669, \
                              21.61736731872703613437482639643069462277, \
                              21.76220974173830183110201170131191655663, \
                              21.90609448143183688097899306297480577485, \
                              22.0490402867972684906486490578908342351, \
                              22.19106530292487557569193393086240056928, \
                              22.33218709789200145155224658654248859863, \
                              22.47242268812972763181726728233480068426, \
                              22.61178856237350107116482179887720435212, \
                              22.75030070429314808911143247007435463568, \
                              22.88797461389019898535743216669742834062, \
                              23.02482532774361176686908350481669578922, \
                              23.16086743817875377107897094013951992889, \
                              23.29611511142881612007056991591400670695, \
                              23.43058210485264430911534250581744912669, \
                              23.56428178326822106049067042363678474115, \
                              23.69722713445669222196499462428128809258, \
                              23.82943078388784484123560074037466926473, \
                              23.96090500871429448226506215187209941985, \
                              24.09166175107828578483607659944185750049, \
                              24.22171263077192876006624179524327743237, \
                              24.35106895728885870286413371602770640852, \
                              24.47974174130269770440759596094085980563, \
                              24.60774170560529058502534730748486421888, \
                              24.73507929553546962640942278460409394305, \
                              24.86176468892705450126529037039701151511, \
                              24.98780780560290158433514050182926229463, \
                              25.11321831644006708870055811559405011212, \
                              25.23800565202952916807855725233848557165, \
                              25.36217901095241434764955105803052194794, \
                              25.48574736769328349128365881159583403571, \
                              25.60871948020974298826286481975975821267, \
                              25.73110389717644977090884369535145447412, \
                              25.85290896492046671466441961895419977548, \
                              25.97414283406389114600188817002274285412, \
                              26.09481346588871741131460938622563226671, \
                              26.21492863843799910251461835808144532308, \
                              26.33449595236654244200346861227048588322, \
                              26.45352283655358479297050162584261891475, \
                              26.57201655348918697331695241490058544311, \
                              26.68998420444539106900228677122542365513, \
                              26.80743273444256315115350431279997534856, \
                              26.92436893702074938634192144279007313002, \
                              27.04079945882532144885899301531481996968, \
                              27.15673080401567010330138314407973778928, \
                              27.27216933850522175671071533349314030512, \
                              27.38712129404059931904444842046619226734, \
                              27.5015927721273236836935924258084318931, \
                              27.61558974780905354225731176234972335349, \
                              27.7291180733069872325213158703418891891, \
                              27.84218348152569918138829195005926719402, \
                              27.95479158943135367217217459578099206443, \
                              28.06694790130792868493048549539313948149, \
                              28.17865781189679108601986330122701795418, \
                              28.28992660942469023626614819841869888214, \
                              28.40075947852497899599227553112094158799, \
                              28.51116150305662806454610751360670020052, \
                              28.62113766882537061488637648812395120886, \
                              28.73069286621109835498093694749057922709, \
                              28.83983189270542661807516764707431009937, \
                              28.94855945536315406496103915669137916186, \
                              29.0568801731711613409334377673415823046, \
                              29.16479857933812188755478470256750080465, \
                              29.27231912350823643173581130305598484097, \
                              29.37944617390204987306787317871788072344, \
                              29.48618401938726481655462124120669710125, \
                              29.59253687148232934120849531844604066321, \
                              29.69850886629544727947065323370919097601, \
                              29.80410406640153686430690657788455184085, \
                              29.90932646265954766612311748559774855157, \
                              30.01417997597243590390627603971350150613, \
                              30.11866845899199411336949914690532327191, \
                              30.22279569777063245220418800300579966128, \
                              30.32656541336211530366113731323828860625, \
                              30.42998126337316800985656774510794595015, \
                              30.53304684346778424965974281784788344284, \
                              30.63576568882598451462550996210811760772, \
                              30.73814127555870008843862044749237208531, \
                              30.84017702208038467423650413525089107265, \
                              30.94187629044088712766532761964247197981, \
                              31.04324238761805344250874245193788967589, \
                              31.14427856677246401343051602820994092765, \
                              31.24498802846565309148851505340619403559, \
                              31.34537392184310108800767554807196747595, \
                              31.44543934578323681652762544616436194505, \
                              31.54518735001363574556992207115040025218, \
                              31.64462093619555173030629520550814059172])
    y,x=0.01*numpy.mgrid[1:101,1:2001]
    c=numpy.zeros_like(x)
    cshort=numpy.zeros_like(x)
    cerf=numpy.zeros_like(x)
    
    num_terms=40
    cwall=1
    c0=0
    coeffm=wmfivetwo/wmcubic*(c0-cwall)
    omega=1.8
    u_bubble=0.05*40
    diffusion=1.0/3.0*(1.0/omega-0.5)

    pe=u_bubble/diffusion
    
    sign=-1
    for i in range(0,num_terms):
        sign=-sign
        cerf=cerf+sign*(scipy.special.erfc((y+2.0*i)/numpy.sqrt(4.0*x/pe))+scipy.special.erfc((2.0*(i+1)-y)/numpy.sqrt(4.0*x/pe)))
    cerf=c0-(c0-cwall)*cerf
    
    for i in range(0,num_terms):
        if i<40:
            cshort=cshort+coeffm[i]*scipy.special.jv(0.25,bessel_zeros[i]*bessel_zeros[i]*0.5*y*y)*numpy.sqrt(y)*numpy.exp(-(bessel_zeros[i]**4)*x*diffusion/u_bubble)
        c=c+coeffm[i]*scipy.special.jv(0.25,bessel_zeros[i]*bessel_zeros[i]*0.5*y*y)*numpy.sqrt(y)*numpy.exp(-(bessel_zeros[i]**4)*x*diffusion/u_bubble)
    c=c+cwall
    cshort=cshort+cwall
    pylab.figure()
    pylab.imshow(c,extent=(0,10,0,0.5))
    pylab.colorbar()
 
    c_levels=numpy.arange(0.2,1.0,0.1)
 
    conc_antibb =numpy.loadtxt(file_dir+"FullProfile/"+"film_antibb0050000.dat")
    conc_inamuro=numpy.loadtxt(file_dir+"FullProfile/"+"film_outflow0050000.dat")
    conc_orest=numpy.loadtxt("/home/shurik/Dropbox/Projects/SharedAlexOrest/Data.txt")
    dims=conc_antibb.shape    
    print "Dims=",dims
    print "Pe=",pe
    
    pylab.figure()    
    pylab.imshow(conc_antibb[0:dims[0]/2,:],extent=(0,20,0,0.5))
    pylab.title("AntiBB simulations")
    pylab.colorbar()    
    pylab.figure()
    pylab.imshow(conc_inamuro[0:dims[0]/2,:],extent=(0,20,0,0.5))
    pylab.title("Inamuro simulations")
    pylab.colorbar()
    pylab.figure()
    orest_dims=conc_orest.shape
    #combined_array=[]
    #for i in range(orest_dims[0]-1):
    #    if conc_orest[i,0]==-100 and conc_orest[i,1]==-100:
    #            continue        
    #    for j in range(i+1,orest_dims[0]):
    #        if conc_orest[i,0]==conc_orest[j,0] and conc_orest[i,1]==conc_orest[j,1]:
    #            conc_orest[j,0]=-100
    #            conc_orest[j,1]=-100
    #    combined_array.append(conc_orest[i])
    #    print i
    #combined_array=numpy.array(combined_array)
    #numpy.savetxt(file_dir+"FullProfile/"+"orest.txt",combined_array)
    
    combined_array=numpy.loadtxt(file_dir+"FullProfile/orest.txt")    

        
#        for (counter2,x2,y2) in enumerate(conc_orest[:,:2]):
#            if x==x2 and y==y2:
#                print counter1,counter2
    triang = tri.Triangulation(combined_array[:,0], combined_array[:,1])
    print triang.triangles.shape
    ctriang=pylab.tricontour(triang,combined_array[:,2],levels=c_levels)
    pylab.clabel(ctriang)
    
    pylab.figure(99,figsize=(20,1))
    c1=pylab.contour(conc_antibb,levels=c_levels,extent=(0.0,20.0,0.0,1.0))
    pylab.clabel(c1,fontsize=9, inline=1)    
    #pylab.figure(figsize=(10,1))
    c2=pylab.contour(conc_inamuro,levels=c_levels,extent=(0.0,20.0,0.0,1.0))
    pylab.clabel(c2,fontsize=9, inline=1)
    
    
    #c_anal=pylab.contour(c,levels=c_levels,linestyles="dashed",extent=(0.0,10.0,0.0,0.5))
    #c_anal_short=pylab.contour(cshort,levels=c_levels,linestyles="dotted",extent=(0.0,10.0,0.0,0.5))
    
    ctriang=pylab.tricontour(triang,combined_array[:,2],levels=c_levels,linestyles="dashdot")
    print cshort
    #c_anal_erf=pylab.contour(cerf,levels=c_levels,linestyles="dotted",extent=(0.0,10.0,0.0,0.5))
    pylab.figure()
    pylab.plot(cshort[dims[0]/2,:])
    pylab.figure()
    xcoor=numpy.linspace(0.0,20.0,dims[1])
    print xcoor    
    print xcoor.shape
    print conc_antibb[dims[0]/2,:].shape
    pylab.plot(xcoor,conc_antibb[dims[0]/2,:])
    pylab.plot(xcoor,conc_inamuro[dims[0]/2,:])

if __name__=="__main__":
    #file_name="../Benchmarks/density0001000.dat"
    file_name="../Benchmarks/conc_initial.dat"    
    file_dir="../Benchmarks/"
    #read_file(file_name)
    
    #read_files(file_dir)
    #show_bessel()    
    #read_film(file_dir)    
    #film_analytical()    
    #compare_film(file_dir)
    #compare_three_films(file_dir)
    #compare_full_profiles(file_dir)    
    film_reconstruction(file_dir)
    pylab.show()
