module emission

implicit none

integer, parameter                  :: Ny=200, dp=8, sp=4
real(dp), parameter                 :: pi=acos(-1.d0), b=6.83089d0, zs=1.6183d-4, &
                                       gg=10.246d0, Q=0.35999d0, kBoltz=8.6173324d-5
logical, save                       :: spectroscopy= .false.

type, public    :: EmissionData

    real(dp)    :: F=5.d0, R=5.d0, gamma=1.d1 
        !Electrostatics: Field, Radius, enhancement factor
    real(dp)    :: W=4.5d0, kT=2.5d-2
        !Materical Characteristics: Work funciton,Temperature
    real(dp)    :: Jem=1.d-200, heat=1.d-200, Gam=1.d10 
        !Calculated results: Current density and Nott. heat
    real(dp)    :: xm=-1.d20, Um=-1.d20, maxbeta=0.d0, minbeta=0.d0
        !Barrier characteristics: xm=x point where barrier is max,
        ! Um=maximum barrier value, dGUm=dG/dW at E=Um, dGEf=dG/dW at Fermi,
    character   :: regime ='F', sharpness = 'B'
        !'f' for field, 'i' for intermediate, 't' for thermal regimes
        !'s' for sharp tip (numerical integration) and 'b' for blunt (KX approx)
    logical     :: full !full calculation if true, else GTF approximation
end type EmissionData


 contains

subroutine cur_dens(this)
!Calculates current density, main module function

    type(EmissionData), intent(inout)      :: this
    
    real(dp)        :: Jf,Jt,n,s,E0,dE,heatf,heatt
    real(dp)        :: nlimf=.7d0, nlimt=2.d0
    
    if (this%R < 0.d0) then
        this%R = 1.d4
    endif
    
    if (this%gamma < 0.d0) then
        this%gamma = 1.d0
    endif
    
    if (this%full) then
        nlimf = 0.7d0
        nlimt = 2.d0
    else
        nlimf = 1.d0
        nlimt = 1.d0
    endif
    
    this%Um = -1.d20
    this%xm = -1.d20
    call gamow_general(this,.true.)

    if (this%kT * this%maxbeta < nlimf) then!field regime
        this%Jem = zs*pi* this%kT * exp(-this%Gam) &
                    / (this%maxbeta * sin(pi * this%maxbeta * this%kT))
        !Murphy-Good version of FN
        this%heat = -zs*(pi**2) * (this%kT**2) * &
                    exp(-this%Gam)*cos(pi* this%maxbeta * this%kT)  &
                    / (this%maxbeta * (sin(pi * this%maxbeta * this%kT )**2))
        this%regime = 'F'
    else if (this%kT * this%minbeta > nlimt) then !thermal regime
        n = 1.d0/(this%kT * this%minbeta)
        s = this%minbeta * this%Um
        Jf = zs * ( this%minbeta **(-2)) * Sigma(1.d0/n) * exp(-s);
        Jt = zs * (this%kT**2) * exp(-n*s) * Sigma(n)
        this%Jem = Jt+Jf/(n**2)
        heatt = (n*s+1.d0)*Sigma(n)*exp(-n*s)
        heatf = (n**3)*DSigma(1.d0/n)*exp(-s)
        this%heat = zs*(heatt-heatf)*(this%kT**3)
        this%regime = 'T'
    else  !intermediate regime
        if (this%full) then
            call J_num_integ(this) !Numerical integration over energies
        else
            call GTFinter(this)
        endif
        this%regime = 'I'
    endif
    
    contains

    pure function Sigma(x) result(Sig)
        double precision, intent(in) :: x
        double precision:: Sig
        Sig = (1.d0+x**2)/(1.d0-x**2)-.3551d0*x**2-0.1059d0*x**4-0.039*x**6!Jensen's sigma
    end function Sigma
    
    pure function DSigma(x) result(ds)
        real(dp), intent(in):: x
        real(dp):: ds
        ds=(-1.d0+4.d0*x**2+x**4)/(1.d0-x**2)-.3551d0*x**2-.3178d0*x**4-0.027066*x**6!Jensen's sigma'
    end function DSigma

end subroutine cur_dens 


subroutine gamow_general(this,full)
!Calculates Gamow exponent in the general case: Choose appropriate approximation

    type(EmissionData), intent(inout)       :: this
    type(EmissionData)                      :: new
    logical, intent(in)                     :: full

    real(dp), parameter                     :: dw = 1.d-2, xlim = 0.08d0
    real(dp)                                :: x, xmax
    
    x = this%W / (this%F * this%R)

    if (x>xlim) then !not x<<1, use numerical integration
        if (x>0.4d0) then !the second root is far away
            xmax = this%W * this%gamma / this%F
        else  !the second root is close to the standard W/F
            xmax = 2.d0 * this%W / this%F
        endif
        call gamow_num(this,xmax)
        this%sharpness='S'
        if (this%Um < 0.d0 .or. (.not. full)) then !barrier lost
            this%maxbeta = 0.d0
        elseif (this%Gam == 1.d20) then
            this%maxbeta = 1.d20
        else
            new = this
            new%W = this%W + dw
            call gamow_num(new,xmax)
            this%maxbeta = abs(new%Gam - this%Gam) / dw
        endif
    else !large radius, KX approximation usable
        call gamow_KX(this)
        this%sharpness='B'
    endif
        
end subroutine gamow_general

subroutine gamow_KX(this) 
    !calculate Gamow parameters using KX approximation
   
    type(EmissionData), intent(inout)   :: this
    real(dp), parameter, dimension(Ny)  :: & !these are the special functions (see Kyritsakis, Xanthakis, PRSA 471:20140811) 
                                           vv = [1.000000000000000e+00, 9.901487508002239e-01, 9.815990523568811e-01,9.735385052587654e-01,9.657944514409587e-01,9.582849898122546e-01,9.509621438046620e-01,9.437939679250140e-01,9.367579341725516e-01,9.298371861358989e-01,9.230186634171222e-01,9.162918241873717e-01,9.096479972736629e-01,9.030804149675133e-01,8.965829830311779e-01,8.901503987642890e-01,8.837783939899669e-01,8.774630622855859e-01,8.712009093011192e-01,8.649890650533810e-01,8.588244941766686e-01,8.527048450122413e-01,8.466279728602170e-01,8.405921698573996e-01,8.345949882169477e-01,8.286351819614325e-01,8.227111386824252e-01,8.168214188018590e-01,8.109647744130352e-01,8.051397093423305e-01,7.993453448084236e-01,7.935805392550657e-01,7.878444892870563e-01,7.821357877770025e-01,7.764539521629619e-01,7.707981061398350e-01,7.651674610258777e-01,7.595613054181773e-01,7.539789715378010e-01,7.484194923946329e-01,7.428826442049090e-01,7.373678107420274e-01,7.318741030001834e-01,7.264014059236763e-01,7.209489040937604e-01,7.155161967568670e-01,7.101029669799113e-01,7.047087080077287e-01,6.993327657426716e-01,6.939749433205432e-01,6.886348537302909e-01,6.833121155653594e-01,6.780063732231738e-01,6.727172947913285e-01,6.674444563060048e-01,6.621875170015923e-01,6.569463113632747e-01,6.517205871223191e-01,6.465099244068481e-01,6.413139630242783e-01,6.361327655966146e-01,6.309656631281146e-01,6.258127673563729e-01,6.206735740822522e-01,6.155481092099058e-01,6.104358084012765e-01,6.053368610270500e-01,6.002506752302079e-01,5.951772367709246e-01,5.901164432886615e-01,5.850679430784230e-01,5.800314714307956e-01,5.750070316723157e-01,5.699944245679107e-01,5.649934599569441e-01,5.600039562421042e-01,5.550257399122128e-01,5.500585920685797e-01,5.451024262683717e-01,5.401571050142602e-01,5.352224820426076e-01,5.302984083461624e-01,5.253846727012907e-01,5.204812022060537e-01,5.155878741643488e-01,5.107045708333732e-01,5.058311791734402e-01,5.009675906125862e-01,4.961136875835652e-01,4.912691037168074e-01,4.864339969693915e-01,4.816082751853399e-01,4.767917683714146e-01,4.719841616151277e-01,4.671856609941015e-01,4.623961033394886e-01,4.576151202714256e-01,4.528429947880911e-01,4.480793569726614e-01,4.433242129835121e-01,4.385775586513992e-01,4.338390551284855e-01,4.291089661228503e-01,4.243867851005731e-01,4.196728691135294e-01,4.149667137737938e-01,4.102686023279039e-01,4.055781952693642e-01,4.008955418859766e-01,3.962206235571218e-01,3.915531020915147e-01,3.868932530612438e-01,3.822407324956104e-01,3.775955424031231e-01,3.729577701211998e-01,3.683271031128482e-01,3.637035566567341e-01,3.590871937782671e-01,3.544778849763434e-01,3.498753504769608e-01,3.452797836412470e-01,3.406911106302504e-01,3.361092595666234e-01,3.315339415893150e-01,3.269653057332668e-01,3.224033030501724e-01,3.178478686121655e-01,3.132989391546318e-01,3.087564443036809e-01,3.042201949462451e-01,2.996902846932734e-01,2.951666561764913e-01,2.906492534459388e-01,2.861380219254216e-01,2.816329083696306e-01,2.771338608228555e-01,2.726408285792251e-01,2.681537621444056e-01,2.636726131986982e-01,2.591973345614714e-01,2.547278801568779e-01,2.502642049807954e-01,2.458062650689477e-01,2.413540174661489e-01,2.369074201966328e-01,2.324664322354155e-01,2.280310134806558e-01,2.236011247269666e-01,2.191767276396464e-01,2.147577847297865e-01,2.103442593302260e-01,2.059361155723163e-01,2.015332821858309e-01,1.971356532159229e-01,1.927432929496386e-01,1.883561686433637e-01,1.839742482442248e-01,1.795975003711649e-01,1.752258942966304e-01,1.708592517296071e-01,1.664976274499445e-01,1.621410476170106e-01,1.577894840762422e-01,1.534429092411601e-01,1.491010897296571e-01,1.447641544395267e-01,1.404321200748444e-01,1.361049612954287e-01,1.317824526245003e-01,1.274646793909447e-01,1.231517008329012e-01,1.188434746806572e-01,1.145396994344216e-01,1.102406427428393e-01,1.059462826904284e-01,1.016563417872831e-01,9.737096326668701e-02,9.309021099940061e-02,8.881381387992358e-02,8.454188344436472e-02,8.027451288980818e-02,7.601136579462732e-02,7.175266600595823e-02,6.749841745663981e-02,6.324828733292226e-02,5.900261251872982e-02,5.476114825994762e-02,5.052390532158921e-02,4.629102797114641e-02,4.206212445720553e-02,3.783758113487127e-02,3.361705982468581e-02,2.940072995229685e-02,2.518850831056418e-02,2.098030715448473e-02,1.677627349067132e-02,1.257611923377084e-02,8.380165347191253e-03,4.187978961787145e-03,0.000000000000000e+00 ], &
                                           tt = [1.000000000000000e+00,1.002030781417862e+00,1.003632786824882e+00,1.005075955514299e+00,1.006417353308358e+00,1.007683957539755e+00,1.008891566532568e+00,1.010050610037040e+00,1.011168454232591e+00,1.012250589592541e+00,1.013301261227889e+00,1.014323859598625e+00,1.015321150800558e+00,1.016295389546883e+00,1.017248511247342e+00,1.018182174987559e+00,1.019097785818562e+00,1.019996584163601e+00,1.020879664647628e+00,1.021747978811433e+00,1.022602411364004e+00,1.023443726618845e+00,1.024272615800069e+00,1.025089688973000e+00,1.025895564002957e+00,1.026690735581673e+00,1.027475690358465e+00,1.028250872863402e+00,1.029016688121873e+00,1.029773532609849e+00,1.030521741770611e+00,1.031261648346635e+00,1.031993546755459e+00,1.032717751209773e+00,1.033434506233224e+00,1.034144066588887e+00,1.034846669618536e+00,1.035542537051420e+00,1.036231877908698e+00,1.036914907534409e+00,1.037591792545195e+00,1.038262712826807e+00,1.038927854335332e+00,1.039587356298646e+00,1.040241387041531e+00,1.040890087747245e+00,1.041533590044200e+00,1.042172029545119e+00,1.042805543516137e+00,1.043434242900508e+00,1.044058243348362e+00,1.044677656087644e+00,1.045292587255446e+00,1.045903138257185e+00,1.046509410981598e+00,1.047111500269847e+00,1.047709490464150e+00,1.048303466845580e+00,1.048893518817778e+00,1.049479730418012e+00,1.050062166197551e+00,1.050640917694896e+00,1.051216043057471e+00,1.051787622582175e+00,1.052355713907147e+00,1.052920395269537e+00,1.053481714557585e+00,1.054039747096986e+00,1.054594545022903e+00,1.055146162284781e+00,1.055694660205253e+00,1.056240095324423e+00,1.056782513399835e+00,1.057321965993506e+00,1.057858503006047e+00,1.058392172758164e+00,1.058923022067008e+00,1.059451098037307e+00,1.059976442316032e+00,1.060499097189592e+00,1.061019104199625e+00,1.061536503973864e+00,1.062051338074067e+00,1.062563642847344e+00,1.063073455326506e+00,1.063580811550365e+00,1.064085746605300e+00,1.064588294664522e+00,1.065088489404378e+00,1.065586370798085e+00,1.066081963100914e+00,1.066575297145695e+00,1.067066405247904e+00,1.067555322913718e+00,1.068042070814106e+00,1.068526679173990e+00,1.069009182922810e+00,1.069489599221223e+00,1.069967961869527e+00,1.070444294259226e+00,1.070918619653641e+00,1.071390969201728e+00,1.071861358501953e+00,1.072329821768765e+00,1.072796371550207e+00,1.073261040898029e+00,1.073723843463047e+00,1.074184807337213e+00,1.074643950999066e+00,1.075101294184407e+00,1.075556863526736e+00,1.076010671316672e+00,1.076462743496344e+00,1.076913097845946e+00,1.077361750009438e+00,1.077808723997548e+00,1.078254036338585e+00,1.078697702306128e+00,1.079139740968321e+00,1.079580174122345e+00,1.080019013476803e+00,1.080456276036589e+00,1.080891978476461e+00,1.081326141457119e+00,1.081758776663107e+00,1.082189899498954e+00,1.082619525426810e+00,1.083047669623051e+00,1.083474347148894e+00,1.083899575009856e+00,1.084323364713017e+00,1.084745730381644e+00,1.085166685890044e+00,1.085586244870236e+00,1.086004420718378e+00,1.086421226600978e+00,1.086836675460887e+00,1.087250780023084e+00,1.087663552800267e+00,1.088075006098254e+00,1.088485152021201e+00,1.088894002476648e+00,1.089301569180396e+00,1.089707863661225e+00,1.090112897265461e+00,1.090516681161385e+00,1.090919226343516e+00,1.091320543636742e+00,1.091720643700332e+00,1.092119537031810e+00,1.092517233970724e+00,1.092913744702279e+00,1.093309079782908e+00,1.093703250103651e+00,1.094096263974417e+00,1.094488130990864e+00,1.094878860611112e+00,1.095268462158765e+00,1.095656944825848e+00,1.096044319655648e+00,1.096430594403337e+00,1.096815777135447e+00,1.097199876540108e+00,1.097582901187470e+00,1.097964862139642e+00,1.098345765703773e+00,1.098725619561639e+00,1.099104431830844e+00,1.099482212946610e+00,1.099858969417514e+00,1.100234708029298e+00,1.100609436704620e+00,1.100983166253060e+00,1.101355900746938e+00,1.101727647595793e+00,1.102098416970978e+00,1.102468214180127e+00,1.102837045403800e+00,1.103204920390876e+00,1.103571844600817e+00,1.103937823680700e+00,1.104302867721525e+00,1.104666980816567e+00,1.105030169345461e+00,1.105392442996972e+00,1.105753804595549e+00,1.106114262768074e+00,1.106473823479837e+00,1.106832491373410e+00,1.107190276156247e+00,1.107547180034785e+00,1.107903212018282e+00,1.108258376274143e+00,1.108612679272688e+00,1.108966127420483e+00,1.109318724982190e+00,1.109670479992468e+00,1.110021395106194e+00,1.110371479424545e+00,1.110720734539592e+00 ], &
                                           ww = [8.000000000000000e-01,7.990593903557062e-01,7.981212218435605e-01,7.971850716193163e-01,7.962506895376646e-01,7.953179118471086e-01,7.943866531707875e-01,7.934567967549424e-01,7.925282750106136e-01,7.916010255446501e-01,7.906749965919331e-01,7.897501209920006e-01,7.888263277835517e-01,7.879036387000090e-01,7.869820002334283e-01,7.860613472701492e-01,7.851416836616130e-01,7.842229778028280e-01,7.833051903036342e-01,7.823883315435012e-01,7.814723156121555e-01,7.805571397772499e-01,7.796427974499822e-01,7.787293322050833e-01,7.778165990366039e-01,7.769046707144984e-01,7.759935170618699e-01,7.750831196172159e-01,7.741734788014050e-01,7.732645222209847e-01,7.723562918756240e-01,7.714487578483593e-01,7.705419490885239e-01,7.696357496513833e-01,7.687302333565893e-01,7.678253773672641e-01,7.669211679764979e-01,7.660175993620985e-01,7.651146669709347e-01,7.642122909578154e-01,7.633105511245557e-01,7.624094327879908e-01,7.615088525828341e-01,7.606088972682540e-01,7.597094886173841e-01,7.588106378635329e-01,7.579123710656784e-01,7.570146661898948e-01,7.561174624569775e-01,7.552208012417980e-01,7.543246762447402e-01,7.534290789679281e-01,7.525340036541388e-01,7.516394470233185e-01,7.507453810173614e-01,7.498517943075966e-01,7.489587146945831e-01,7.480661470849971e-01,7.471740538712757e-01,7.462824091535802e-01,7.453912867492993e-01,7.445005821810872e-01,7.436103774063627e-01,7.427206035994978e-01,7.418313187646799e-01,7.409424355359369e-01,7.400540488502052e-01,7.391660610822821e-01,7.382785145323230e-01,7.373914287447291e-01,7.365047604128426e-01,7.356184861552344e-01,7.347326478010596e-01,7.338472359139018e-01,7.329622421200244e-01,7.320776590403786e-01,7.311934802272851e-01,7.303096867745417e-01,7.294262920302131e-01,7.285432956473210e-01,7.276606941730114e-01,7.267784825772383e-01,7.258966391623030e-01,7.250151764682804e-01,7.241340936738493e-01,7.232533904789240e-01,7.223730670724802e-01,7.214931241023815e-01,7.206135592443398e-01,7.197343054481833e-01,7.188554297854226e-01,7.179769346016311e-01,7.170988015262763e-01,7.162209738905005e-01,7.153435294366657e-01,7.144664498441585e-01,7.135896626254823e-01,7.127132644162992e-01,7.118371812429732e-01,7.109614366772887e-01,7.100860513209014e-01,7.092109574472614e-01,7.083362451607804e-01,7.074618015503424e-01,7.065877407427119e-01,7.057139492352485e-01,7.048405211933729e-01,7.039673856690785e-01,7.030945743065227e-01,7.022221004806966e-01,7.013498921102026e-01,7.004780402547411e-01,6.996064704755764e-01,6.987352002768673e-01,6.978642698311666e-01,6.969936118137134e-01,6.961232464350917e-01,6.952532065889611e-01,6.943834731280344e-01,6.935139859278413e-01,6.926448122914752e-01,6.917759472301522e-01,6.909073860082978e-01,6.900390646592811e-01,6.891710377861378e-01,6.883033061043673e-01,6.874358657452273e-01,6.865687130507158e-01,6.857018421806050e-01,6.848352145252686e-01,6.839688686162471e-01,6.831028015145367e-01,6.822370104572679e-01,6.813714928513894e-01,6.805062462675940e-01,6.796412684344860e-01,6.787765572329649e-01,6.779121106908332e-01,6.770479269776047e-01,6.761840043995144e-01,6.753203413947115e-01,6.744569365286324e-01,6.735937884895586e-01,6.727308960843218e-01,6.718682582341828e-01,6.710058739708529e-01,6.701437424326648e-01,6.692818628608863e-01,6.684202345961574e-01,6.675588570750687e-01,6.666977298268534e-01,6.658368524702023e-01,6.649762145174657e-01,6.641157955198331e-01,6.632556228842422e-01,6.623956965898202e-01,6.615360166902944e-01,6.606765833114356e-01,6.598173966485905e-01,6.589584148830649e-01,6.580996621472397e-01,6.572411544304840e-01,6.563828922121562e-01,6.555248760295980e-01,6.546670475669887e-01,6.538094516821134e-01,6.529521012357424e-01,6.520949970044247e-01,6.512380822956181e-01,6.503813892668515e-01,6.495249426319587e-01,6.486687379272718e-01,6.478126959530650e-01,6.469569010900136e-01,6.461013544668252e-01,6.452459833377455e-01,6.443908361595955e-01,6.435359385573338e-01,6.426812191722087e-01,6.418267171808800e-01,6.409724665326726e-01,6.401183764547567e-01,6.392645186054525e-01,6.384109008745196e-01,6.375574330096243e-01,6.367042195871654e-01,6.358511957912952e-01,6.349983711019559e-01,6.341457940200761e-01,6.332933562857685e-01,6.324411770853480e-01,6.315891635727651e-01,6.307373714298024e-01,6.298857823014392e-01,6.290343763414623e-01,6.281832026705599e-01,6.273321823616571e-01,6.264814156974592e-01,6.256307808321395e-01,6.247804132000000e-01 ], &
                                           psi = [1.333333333333333e+00,1.333018024233115e+00,1.332701031413377e+00,1.332382777718021e+00,1.332063462481988e+00,1.331743219800990e+00,1.331422162060289e+00,1.331100359883384e+00,1.330777880951771e+00,1.330454780072042e+00,1.330131104829534e+00,1.329806888425339e+00,1.329482155570348e+00,1.329156960054474e+00,1.328831322965915e+00,1.328505256938984e+00,1.328178797347004e+00,1.327851963327141e+00,1.327524768508933e+00,1.327197243420925e+00,1.326869379644892e+00,1.326541198974261e+00,1.326212720264901e+00,1.325883980612838e+00,1.325554942466070e+00,1.325225651926277e+00,1.324896114017230e+00,1.324566337390244e+00,1.324236337281653e+00,1.323906099425770e+00,1.323575653936342e+00,1.323245001921238e+00,1.322914167210422e+00,1.322583115490618e+00,1.322251887322098e+00,1.321920484386648e+00,1.321588911526883e+00,1.321257176307364e+00,1.320925286409709e+00,1.320593218868366e+00,1.320261014491362e+00,1.319928675858194e+00,1.319596177472917e+00,1.319263562181622e+00,1.318930805941319e+00,1.318597920571705e+00,1.318264923692205e+00,1.317931813205207e+00,1.317598570995208e+00,1.317265220239226e+00,1.316931764557958e+00,1.316598206480150e+00,1.316264549464544e+00,1.315930797806357e+00,1.315596945492719e+00,1.315262993168529e+00,1.314928957364734e+00,1.314594845117375e+00,1.314260645869617e+00,1.313926353694029e+00,1.313592003542662e+00,1.313257556921115e+00,1.312923051920728e+00,1.312588464442063e+00,1.312253822540088e+00,1.311919094199020e+00,1.311584322552062e+00,1.311249471102565e+00,1.310914561069798e+00,1.310579604205516e+00,1.310244586153102e+00,1.309909500694144e+00,1.309574368563633e+00,1.309239189180093e+00,1.308903962321453e+00,1.308568688100718e+00,1.308233366943394e+00,1.307897994020324e+00,1.307562577848092e+00,1.307227121189610e+00,1.306891625453616e+00,1.306556091330614e+00,1.306220512499999e+00,1.305884896854858e+00,1.305549246667971e+00,1.305213564377872e+00,1.304877852577373e+00,1.304542114002722e+00,1.304206350100934e+00,1.303870535200850e+00,1.303534699659165e+00,1.303198846741950e+00,1.302862971001015e+00,1.302527050894820e+00,1.302191121124099e+00,1.301855176135410e+00,1.301519187617303e+00,1.301183198107724e+00,1.300847178591757e+00,1.300511140910904e+00,1.300175095640288e+00,1.299839016245855e+00,1.299502942413006e+00,1.299166828556455e+00,1.298830724405211e+00,1.298494584040807e+00,1.298158448766632e+00,1.297822290469760e+00,1.297486124129521e+00,1.297149957023230e+00,1.296813760453559e+00,1.296477574324301e+00,1.296141368885803e+00,1.295805153046724e+00,1.295468945241037e+00,1.295132718620175e+00,1.294796483164058e+00,1.294460254178512e+00,1.294124025028758e+00,1.293787771784801e+00,1.293451524188628e+00,1.293115281504094e+00,1.292779043080147e+00,1.292442783268936e+00,1.292106526392454e+00,1.291770274037186e+00,1.291434025837305e+00,1.291097781496767e+00,1.290761539780119e+00,1.290425285605418e+00,1.290089036421810e+00,1.289752792165304e+00,1.289416552829235e+00,1.289080318461995e+00,1.288744089164848e+00,1.288407865089826e+00,1.288071646437702e+00,1.287735433456070e+00,1.287399226437465e+00,1.287063025717571e+00,1.286726831673501e+00,1.286390644722119e+00,1.286054465318477e+00,1.285718293954241e+00,1.285382131156234e+00,1.285045977485010e+00,1.284709833533460e+00,1.284373699925540e+00,1.284037577314941e+00,1.283701466383906e+00,1.283365367842030e+00,1.283029282425125e+00,1.282693206586771e+00,1.282357132559316e+00,1.282021072787891e+00,1.281685028116350e+00,1.281348999409467e+00,1.281012987552014e+00,1.280676993447894e+00,1.280341000234313e+00,1.280005018911602e+00,1.279669057029062e+00,1.279333115571353e+00,1.278997195538133e+00,1.278661273046611e+00,1.278325367809626e+00,1.277989486007780e+00,1.277653628705530e+00,1.277317772671696e+00,1.276981932216333e+00,1.276646118494417e+00,1.276310330321828e+00,1.275974534925357e+00,1.275638768649462e+00,1.275303032648981e+00,1.274967296856635e+00,1.274631582411728e+00,1.274295900796550e+00,1.273960222510329e+00,1.273624564747406e+00,1.273288942483700e+00,1.272953317999347e+00,1.272617722187058e+00,1.272282158991539e+00,1.271946590898323e+00,1.271611062678912e+00,1.271275547550641e+00,1.270940050105779e+00,1.270604591406306e+00,1.270269126311487e+00,1.269933705718180e+00,1.269598290985304e+00,1.269262906179969e+00,1.268927544102622e+00,1.268592196921197e+00,1.268256885909111e+00,1.267921578286122e+00,1.267586316905675e+00,1.267251050867749e+00,1.266915967611353e+00 ]
    real(dp)                            :: yf, t, ps, v, omeg

    yf = 1.44d0 * this%F / (this%W ** 2)
    if (yf>1.d0 .or. yf<0.d0) then
        print *, 'Error: yf out of interpolation bounds'
        return
    endif
    
    t = lininterp(tt,0.d0,1.d0,yf)
    ps = lininterp(psi,0.d0,1.d0,yf)
    v = lininterp(vv,0.d0,1.d0,yf)
    omeg = lininterp(ww,0.d0,1.d0,yf)
    
    this%maxbeta = gg*(sqrt(this%W) / this%F)*(t+ps*(this%W / (this%F * this%R)))
    this%Gam = (2*gg/3)*((this%W **1.5d0) / this%F) *(v+omeg*(this%W / (this%F * this%R)))
    this%Um = this%W - 2.d0*sqrt(this%F * Q)-0.75d0 * Q / this%R
    this%minbeta = 16.093d0/sqrt(this%F **1.5d0 / sqrt(Q) - 4* this%F / this%R)
    
end subroutine gamow_KX

subroutine gamow_num(this,xmax)

    type(EmissionData), intent(inout)   :: this
    real(dp), intent(in)                :: xmax
    
    integer, parameter                  :: maxint = 256
    real(dp), parameter                 :: AtolRoot=1.d-6, RtolRoot=1.d-6, &
                                           AtolInt=1.d-5, RtolInt=1.d-5
    
    real(dp), dimension(maxint)         :: ALIST, BLIST, RLIST, ELIST
    real(dp)                            :: ABSERR
    integer                             :: IFLAG, NEVAL, IER, IORD(maxint), LAST
    real(dp)                            :: G, x = 1.d-5 , x1(2), x2(2)
    
    if (this%Um == -1.d20) then ! Um is not initialized
        this%Um = -local_min(x, xmax, 1.d-8, 1.d-8, neg_bar, this%xm)
    endif
    if (this%Um<0.d0) then
        this%Gam=0.d0
        return
    endif
    if (abs(this%Um - (this%W - this%F * this%R))<.2d0) then ! almost without maximum
        this%minbeta = 1.d20
    else
        this%minbeta = 22.761d0 / sqrt(abs(diff2(bar,this%xm)))
    endif
    
    x1 = [0.01d0,this%xm]
    x2 = [this%xm,xmax]
    call dfzero(bar,x1(1),x1(2),x1(1),RtolRoot,AtolRoot,IFLAG)
    call dfzero(bar,x2(1),x2(2),x2(1),RtolRoot,AtolRoot,IFLAG)
    call dqage(sqrt_bar,x1(2),x2(1),AtolInt,RtolInt,2,1000,G,ABSERR, &
    NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
    this%Gam = gg * G
    
    contains
    pure function bar(x) result(V)!sphere barrier model
        real(dp), intent(in)    :: x
        real(dp)                ::V
        V = this%W - (this%F * this%R * x*(this%gamma - 1.d0) + this%F * x**2) &
            / (this%gamma * x + this%R * (this%gamma - 1.d0)) &
            - Q / (x + ( 0.5d0 * (x**2)) / this%R)
    end function bar
    pure function neg_bar(x) result(nv)
        real(dp), intent(in)    ::x
        real(dp)                ::nv
        nv = - bar(x)
    end function neg_bar
    pure function sqrt_bar(x) result(rv)
        real(dp), intent(in)    ::x
        real(dp)                ::rv
        rv=sqrt(bar(x))
    end function sqrt_bar
    
end subroutine gamow_num

subroutine GTFinter(this)
!calculates estimation of current according to GTF theory for intermediate regime
    type(EmissionData), intent(inout)   :: this
    
    real(sp)            :: C_FN, Bq, B_FN, UmaxokT, polybeta(3), work(8)
    complex(sp)         :: rt(2)
    integer             :: ierr
    real(dp)            :: zmax, Em, Gm, s, heatt, heatf
    

    C_FN = real(this%maxbeta * this%Um, sp)
    Bq = real(this%minbeta * this%Um, sp)
    B_FN = real(this%Gam, sp)
    UmaxokT = real(this%Um / this%kT, sp)
    
    polybeta = [3.*(C_FN+Bq-2.*B_FN), 2.*(3.*B_FN-Bq-2.*C_FN), C_FN-UmaxokT]
    
    call RPQR79(2,polybeta,rt,ierr,work)
    zmax=minval(real(rt,dp));
    if (zmax<0) then
        zmax=maxval(real(rt,dp))
    endif
    Em = this%Um * zmax
    Gm = -(C_FN+Bq-2.*B_FN) * zmax**3 - (3.*B_FN-Bq-2.*C_FN) * zmax**2 &
            -C_FN * zmax + B_FN
    s = Em / this%kT + Gm
    
    heatt = (s+1.d0)*exp(-s)
    heatf = exp(-s)
    
    this%heat = zs*(heatt+heatf)*(this%kT**3)
    this%Jem = zs * (this%kT**2) * exp(-s) * (s + 1.d0)

end subroutine GTFinter


subroutine J_num_integ(this)
!numerical integration over energies to obtain total current according
!to landauer formula
    type(EmissionData), intent(inout)   :: this
    type(EmissionData)                  :: new
    
    real(dp), parameter                 :: cutoff=1.d-4 
            !cutoff of the exponentially decreasing
    integer, parameter                  ::Nvals=256, fidout=1953
            !no of intervals for integration and fidout for spectroscopy
            
    real(dp)                            :: Gj, G(4*Nvals), Ej, dE, Ea, Eb
    real(dp)                            :: insum, Jsum,Jcur,Umax,Emax,Emin,E0
    real(dp)                            :: Gup(2*Nvals),Gdown(2*Nvals),fj,fmax
    real(dp)                            :: integrand,outsum
    integer                             :: j, i, k

    if (spectroscopy) then
        open(fidout,file='spectra.csv')
    endif
    this%Um = -1.d20
    this%xm = -1.d20
    
    call gamow_general(this,.true.)
    
    if (this%Gam == 0.d0) then!if no barrier for fermi electrons
        this%Jem = 1.d20
        return
    endif

    if (this%kT * this%maxbeta < 1.05d0) then!field regime, max is close to Ef
        Emax = abs(log(cutoff)) / (1.d0 / this%kT - .7 * this%maxbeta)
        Emin = -abs(log(cutoff)) / this%maxbeta
        E0 = 0.d0
    else if (this%kT * this%minbeta > .95d0) then!thermal regime, max is close to Um
        Emax = this%Um + abs(log(cutoff)) * this%kT
        Emin = this%Um - abs(log(cutoff)) / (this%minbeta - 0.7 / this%kT)
        E0 = this%Um
    else! integmediate regime, max is in between ,find it
        Ea = 0.d0
        Eb = this%Um
        do j=1,100!bisection to find where n=1
            Ej=(Eb+Ea)*.5d0
            Umax = this%Um - Ej
            new = this
            new%W = this%W - Ej
            call gamow_general(new,.true.)
            if (this%kT * new%maxbeta < 0.98d0) then
                Eb = Ej
            else if (this%kT * new%maxbeta > 1.02d0) then
                Ea = Ej
            else
                exit
            endif
        enddo
        Emax = this%Um
        Emin = 0.d0
        E0 = Ej
    endif
    
    dE = (Emax-Emin) / (Nvals-1)
    Jsum = 0.d0
    insum = 0.d0
    outsum = 0.d0
    fmax = 0.d0
    
    do j = 1, 2*Nvals !integrate over high energies
        Ej = E0 + (j-1) * dE
        Umax = this%Um - Ej
        if (Umax > 0.d0) then
            new = this
            new%W = this%W - Ej
            call gamow_general(new,.false.)
            Gj = new%Gam
        else
            Gj = this%minbeta * Umax
        endif
        if ( .not. isnan(this%Gam) ) then
            Gup(j)=Gj
        endif
        fj=lFD(Ej,this%kT)/(1.d0+exp(Gup(j)))
        if (fj>fmax) then
            fmax=fj
        endif
        if (abs(fj/fmax)<cutoff)  then
            exit
        endif
        Jsum = Jsum + fj
    enddo
    
    do i=1,2*Nvals!integrate over low energies
        Ej = E0 - i * dE
        Umax = this%Um - Ej
        if (Umax>0.d0) then
            new = this
            new%W = this%W - Ej
            call gamow_general(new,.false.)
            Gj = new%Gam
        else
            Gj = this%minbeta * Umax
        endif
        if ( .not. isnan(this%Gam) ) then
            Gdown(i)=Gj
        endif
        fj=lFD(Ej,this%kT)/(1.d0+exp(Gdown(i)))
        if (fj>fmax) then
            fmax=fj
        endif
        if (abs(fj/fmax)<cutoff)  then
            exit
        endif
        Jsum = Jsum + fj
    enddo
    this%Jem = Jsum * zs * this%kT * dE
    
    G(1:i+j) = [Gdown(i:1:-1),Gup(1:j)]!bring G vectors together and in correct order
    do k=1,i+j!integrate for nottingham effect
        integrand = Ej * (insum + .5d0*dE / (1.d0 + exp(G(k)))) &
                    / (1.d0 + exp(Ej/this%kT))
        insum = insum + dE / (1.d0 + exp(G(k)))
        outsum = outsum + integrand
        if (spectroscopy) then
            write(fidout,*) Ej,',',integrand,',',integrand/Ej
        endif
        Ej = Ej + dE
    enddo

    this%heat = zs * outsum * dE
    
    if (spectroscopy) then
        close(fidout)
    endif

    contains
    pure function lFD(E,kT) result(L)
        real(dp), intent(in)    ::E,kT
        real(dp)                :: L
        if (E>10.d0*kT) then
            L=exp(-E/kT)
        else
            L=log(1.d0+exp(-E/kT))
        endif
    end function lFD
end subroutine J_num_integ

subroutine print_data(this, filenum)
    !print the state of the object nicely
    type(EmissionData), intent(in)      :: this
    integer, intent(in), optional       :: filenum
    integer                             :: fid
    
    if (.not. present(filenum)) then
        fid = 6
    else
        fid=filenum
    endif
    
    write (fid,'(A10,ES12.4,/A10,ES12.4,/A10,ES12.4)') &
            'F = ', this%F, 'R = ', this%R, 'gamma = ', this%gamma
    write (fid,'(A10,ES12.4,/A10,ES12.4)') 'W = ', this%W, 'kT = ', this%kT
    write (fid,'(A10,ES12.4,/A10,ES12.4,/A10,ES12.4)') &
            'Jem = ', this%Jem, 'heat = ', this%heat, 'Gamow = ', this%Gam
    write (fid,'(A10,ES12.4,/A10,ES12.4,/A10,ES12.4,/A10,ES12.4)') &
            'xm = ', this%xm, 'Um = ', this%Um, &
            'maxbeta = ', this%maxbeta, 'minbeta = ', this%minbeta
    write (fid,'(A10,A12,/A10,A12)') 'Regime:  ', this%regime, &
                                     'Sharpness:', this%sharpness
end subroutine print_data

pure function linspace(a,b,N) result(x)
    real(dp),intent(in)     ::a,b
    integer,intent(in)      ::N
    real(dp)                :: dx,x(N)
    integer                 :: i
    
    dx=(b-a)/(N-1)
    do i=1,N
        x(i)=a+(i-1)*dx
    enddo
end function linspace

function diff2(f,x) result(y)!second derivative of a function f at point x
    real(dp), intent(in)    ::x
    real(dp), external      ::f
    real(dp)                ::y
    real(dp),parameter      ::dx=1.d-2

    y=(f(x+dx)+f(x-dx)-2.d0*f(x))/(dx*dx)
end function diff2

pure function lininterp(yi,a,b,x) result(y)
!simple linear interpolation function
!appropriate for uniform linspace
    real(dp), intent(in)    :: a,b,x, yi(:) 
    !yi interpolation array, a,b are x interval limits and x the requested point
    integer                 :: Nnear,N
    real(dp)                :: y,dx,dy,xnear
    
    if (x<a .or. x>b) then
        y=1.d200
        return
    endif
    N=size(yi)
    Nnear=nint((x-a)*(N-1)/(b-a))+1
    dx=(b-a)/dfloat(N-1)
    xnear=a+(Nnear-1)*dx
    if (x>xnear) then
        dy=yi(Nnear+1)-yi(Nnear)
    else
        dy=yi(Nnear)-yi(Nnear-1)
    endif
    y=yi(Nnear)+(dy/dx)*(x-xnear)
end function

function local_min ( a, b, eps, t, f, x )

!*****************************************************************************80
!! LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
!
!  Discussion:
!    The method used is a combination of golden section search and
!    successive parabolic interpolation.  Convergence is never much slower
!    than that for a Fibonacci search.  If F has a continuous second
!    derivative which is positive at the minimum (which is not at A or
!    B), then convergence is superlinear, and usually of the order of
!    about 1.324....
!
!    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
!    F is never evaluated at two points closer than TOL.  
!
!    If F is a unimodal function and the computed values of F are always
!    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
!    LOCAL_MIN approximates the abscissa of the global minimum of F on the 
!    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.  
!
!    If F is not unimodal, then LOCAL_MIN may approximate a local, but 
!    perhaps non-global, minimum to the same accuracy.
!
!    Thanks to Jonathan Eggleston for pointing out a correction to the 
!    golden section step, 01 July 2013.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!    01 July 2013
!
!  Author:
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!    Richard Brent,
!    Algorithms for Minimization Without Derivatives,
!    Dover, 2002,
!    ISBN: 0-486-41998-3,
!    LC: QA402.5.B74.
!
!  Parameters:
!    Input, real ( kind = 8 ) A, B, the endpoints of the interval.
!
!    Input, real ( kind = 8 ) EPS, a positive relative error tolerance.
!    EPS should be no smaller than twice the relative machine precision,
!    and preferably not much less than the square root of the relative
!    machine precision.
!
!    Input, real ( kind = 8 ) T, a positive absolute error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose local minimum is being sought.
!
!    Output, real ( kind = 8 ) X, the estimated value of an abscissa
!    for which F attains a local minimum value in [A,B].
!
!    Output, real ( kind = 8 ) LOCAL_MIN, the value F(X).
!*****************************************************************************
  implicit none

  real(dp),intent(in)   :: a,b,eps,t
  real(dp),external     :: f
  real(dp),intent(out)  :: x
  real(dp)              :: c,d,e,fu,fv,fw,fx,local_min,m,p,q,r,sa,sb,t2,tol,u,v,w
!
!  C is the square of the inverse of the golden ratio.
!
  c = 0.5D+00 * ( 3.0D+00 - sqrt ( 5.0D+00 ) )

  sa = a
  sb = b
  x = sa + c * ( b - a )
  w = x
  v = w
  e = 0.0D+00
  fx = f ( x )
  fw = fx
  fv = fw

  do
    m = 0.5D+00 * ( sa + sb ) 
    tol = eps * abs ( x ) + t
    t2 = 2.0D+00 * tol
!
!  Check the stopping criterion.
!
    if ( abs ( x - m ) <= t2 - 0.5D+00 * ( sb - sa ) ) then
      exit
    end if
!
!  Fit a parabola.
!
    r = 0.0D+00
    q = r
    p = q

    if ( tol < abs ( e ) ) then

      r = ( x - w ) * ( fx - fv )
      q = ( x - v ) * ( fx - fw )
      p = ( x - v ) * q - ( x - w ) * r
      q = 2.0D+00 * ( q - r )

      if ( 0.0D+00 < q ) then
        p = - p
      end if

      q = abs ( q )

      r = e
      e = d

    end if

    if ( abs ( p ) < abs ( 0.5D+00 * q * r ) .and. &
         q * ( sa - x ) < p .and. &
         p < q * ( sb - x ) ) then
!
!  Take the parabolic interpolation step.
!
      d = p / q
      u = x + d
!
!  F must not be evaluated too close to A or B.
!
      if ( ( u - sa ) < t2 .or. ( sb - u ) < t2 ) then

        if ( x < m ) then
          d = tol
        else
          d = - tol
        end if

      end if
!
!  A golden-section step.
!
    else

      if ( x < m ) then
        e = sb - x
      else
        e = sa - x
      end if

      d = c * e

    end if
!
!  F must not be evaluated too close to X.
!
    if ( tol <= abs ( d ) ) then
      u = x + d
    else if ( 0.0D+00 < d ) then
      u = x + tol
    else
      u = x - tol
    end if

    fu = f ( u )
!
!  Update A, B, V, W, and X.
!
    if ( fu <= fx ) then

      if ( u < x ) then
        sb = x
      else
        sa = x
      end if

      v = w
      fv = fw
      w = x
      fw = fx
      x = u
      fx = fu

    else

      if ( u < x ) then
        sa = u
      else
        sb = u
      end if

      if ( fu <= fw .or. w == x ) then
        v = w
        fv = fw
        w = u
        fw = fu
      else if ( fu <= fv .or. v == x .or. v == w ) then
        v = u
        fv = fu
      end if

    end if

  end do

  local_min = fx

  return
end function local_min

end module emission



!**********************************************************************
!* PROCEDURES FOR THE LEAST-SQUARES SOLUTION OF M NONLINEAR EQUATIONS *
!* IN N VARIABLES USING THE Levenberg-Marquardt ALGORITHM.            *
!* ------------------------------------------------------------------ *
!*  REFERENCE                                                         *
!*  From F77 program By:                                              *
!*  ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. MARCH 1980.         *
!*  BURTON S. GARBOW, KENNETH E. HILLSTROM, JORGE J. MORE             *
!*                                                                    *
!*                               F90 Release By J-P Moreau, Paris     *
!*                                  (See Demo program tlm.f90).       *
!*                                      (www.jpmoreau.fr)             *
!**********************************************************************
MODULE Levenberg_Marquardt
! 25 October 2001:
!    Changed INTENT of iflag in several places to IN OUT.
!    Changed INTENT of fvec to IN OUT in user routine FCN.
!    Replaced several DO loops with array operations.
! amiller @ bigpond.net.au

use emission, only: dp

IMPLICIT NONE
INTEGER, PARAMETER :: Nmaxval = 200

PUBLIC :: lmdif1, lmdif, lmder1, lmder, enorm, nlinfit

 CONTAINS

function nlinfit(fun,xdata,ydata,p0) result(var)
!Added to module By Andreas Kyritsakis 31.03.2016
!Simple function to data to non linear function

    real(dp), intent(in)        :: xdata(:) !input x data
    real(dp), intent(in)        :: ydata(:) !input y data
    real(dp), intent(inout)     :: p0(:)    !in: initial guess, out:result
    interface                               !user provided function to be fit
        pure function fun(x,p) result(y)
            implicit none
            integer,parameter :: dp = selected_real_kind(12,60)
            real(dp), intent(in) :: x
            real(dp), intent(in) :: p(:)
            real(dp)             :: y
        end function fun
    end interface

    real(dp)                    :: var,fvec(size(xdata))
    real(dp),parameter          :: tol=1.d-10 !change tolerance if needed
    integer                     :: i,m,n,info,iwa(size(p0))
    
    m=size(xdata)
    n=size(p0)

    call lmdif1(fcn,m,n,p0,fvec,tol,info,iwa)
    var=sum(sqrt(fvec)/abs(ydata))/m
    
    contains

    subroutine fcn(m, n, p, fvec, iflag)
        integer, intent(in)     :: m, n
        real(dp), intent(in)    :: p(:)
        real(dp), intent(inout) :: fvec(:)
        integer, intent(inout)  :: iflag

        do i=1,m
            fvec(i)=(fun(xdata(i),p)-ydata(i))**2
        enddo
    end subroutine  
end function nlinfit

SUBROUTINE lmdif1(fcn, m, n, x, fvec, tol, info, iwa)
!USE common_refnum 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-11  Time: 00:51:44

! N.B. Arguments WA & LWA have been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(OUT)     :: fvec(:)
REAL (dp), INTENT(IN)      :: tol
INTEGER, INTENT(OUT)       :: info
INTEGER, INTENT(OUT)       :: iwa(:)

! EXTERNAL fcn

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn
END INTERFACE

!  **********

!  subroutine lmdif1

!  The purpose of lmdif1 is to minimize the sum of the squares of m nonlinear
!  functions in n variables by a modification of the Levenberg-Marquardt
!  algorithm.  This is done by using the more general least-squares
!  solver lmdif.  The user must provide a subroutine which calculates the
!  functions.  The jacobian is then calculated by a forward-difference
!  approximation.

!  the subroutine statement is

!    subroutine lmdif1(fcn, m, n, x, fvec, tol, info, iwa)

!  where

!    fcn is the name of the user-supplied subroutine which calculates
!      the functions.  fcn must be declared in an external statement in the
!      user calling program, and should be written as follows.

!      subroutine fcn(m, n, x, fvec, iflag)
!      integer m, n, iflag
!      REAL (dp) x(n), fvec(m)
!      ----------
!      calculate the functions at x and return this vector in fvec.
!      ----------
!      return
!      end

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of lmdif1.
!      In this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number of functions.

!    n is a positive integer input variable set to the number of variables.
!      n must not exceed m.

!    x is an array of length n.  On input x must contain an initial estimate
!      of the solution vector.  On output x contains the final estimate of
!      the solution vector.

!    fvec is an output array of length m which contains
!      the functions evaluated at the output x.

!    tol is a nonnegative input variable.  Termination occurs when the
!      algorithm estimates either that the relative error in the sum of
!      squares is at most tol or that the relative error between x and the
!      solution is at most tol.

!    info is an integer output variable.  If the user has terminated execution,
!      info is set to the (negative) value of iflag.  See description of fcn.
!      Otherwise, info is set as follows.

!      info = 0  improper input parameters.

!      info = 1  algorithm estimates that the relative error
!                in the sum of squares is at most tol.

!      info = 2  algorithm estimates that the relative error
!                between x and the solution is at most tol.

!      info = 3  conditions for info = 1 and info = 2 both hold.

!      info = 4  fvec is orthogonal to the columns of the
!                jacobian to machine precision.

!      info = 5  number of calls to fcn has reached or exceeded 200*(n+1).

!      info = 6  tol is too small. no further reduction in
!                the sum of squares is possible.

!      info = 7  tol is too small.  No further improvement in
!                the approximate solution x is possible.

!    iwa is an integer work array of length n.

!    wa is a work array of length lwa.

!    lwa is a positive integer input variable not less than m*n+5*n+m.

!  subprograms called

!    user-supplied ...... fcn

!    minpack-supplied ... lmdif

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  **********
INTEGER   :: maxfev, mode, nprint
INTEGER, SAVE :: nprob, nfev, njev
REAL (dp) :: epsfcn, ftol, gtol, xtol, fjac(m,n)
REAL (dp), PARAMETER :: factor = 100._dp, zero = 0.0_dp

info = 0

!     check the input parameters for errors.

IF (n <= 0 .OR. m < n .OR. tol < zero) GO TO 10

!     call lmdif.

maxfev = Nmaxval*(n + 1)
ftol = tol
xtol = tol
gtol = zero
epsfcn = zero
mode = 1
nprint = 0
CALL lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,   &
           mode, factor, nprint, info, nfev, fjac, iwa)
IF (info == 8) info = 4

10 RETURN

!     last card of subroutine lmdif1.

END SUBROUTINE lmdif1



SUBROUTINE lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,  &
                 mode, factor, nprint, info, nfev, fjac, ipvt)
 
! N.B. Arguments LDFJAC, DIAG, QTF, WA1, WA2, WA3 & WA4 have been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(OUT)     :: fvec(:)
REAL (dp), INTENT(IN)      :: ftol
REAL (dp), INTENT(IN)      :: xtol
REAL (dp), INTENT(IN OUT)  :: gtol
INTEGER, INTENT(IN OUT)    :: maxfev
REAL (dp), INTENT(IN OUT)  :: epsfcn
INTEGER, INTENT(IN)        :: mode
REAL (dp), INTENT(IN)      :: factor
INTEGER, INTENT(IN)        :: nprint
INTEGER, INTENT(OUT)       :: info
INTEGER, INTENT(OUT)       :: nfev
REAL (dp), INTENT(OUT)     :: fjac(:,:)    ! fjac(ldfjac,n)
INTEGER, INTENT(OUT)       :: ipvt(:)

! EXTERNAL fcn

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn
END INTERFACE

!  **********

!  subroutine lmdif

!  The purpose of lmdif is to minimize the sum of the squares of m nonlinear
!  functions in n variables by a modification of the Levenberg-Marquardt
!  algorithm.  The user must provide a subroutine which calculates the
!  functions.  The jacobian is then calculated by a forward-difference
!  approximation.

!  the subroutine statement is

!    subroutine lmdif(fcn, m, n, x, fvec, ftol, xtol, gtol, maxfev, epsfcn,
!                     diag, mode, factor, nprint, info, nfev, fjac,
!                     ldfjac, ipvt, qtf, wa1, wa2, wa3, wa4)

! N.B. 7 of these arguments have been removed in this version.

!  where

!    fcn is the name of the user-supplied subroutine which calculates the
!      functions.  fcn must be declared in an external statement in the user
!      calling program, and should be written as follows.

!      subroutine fcn(m, n, x, fvec, iflag)
!      integer m, n, iflag
!      REAL (dp) x(:), fvec(m)
!      ----------
!      calculate the functions at x and return this vector in fvec.
!      ----------
!      return
!      end

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of lmdif.
!      in this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number of functions.

!    n is a positive integer input variable set to the number of variables.
!      n must not exceed m.

!    x is an array of length n.  On input x must contain an initial estimate
!      of the solution vector.  On output x contains the final estimate of the
!      solution vector.

!    fvec is an output array of length m which contains
!      the functions evaluated at the output x.

!    ftol is a nonnegative input variable.  Termination occurs when both the
!      actual and predicted relative reductions in the sum of squares are at
!      most ftol.  Therefore, ftol measures the relative error desired
!      in the sum of squares.

!    xtol is a nonnegative input variable.  Termination occurs when the
!      relative error between two consecutive iterates is at most xtol.
!      Therefore, xtol measures the relative error desired in the approximate
!      solution.

!    gtol is a nonnegative input variable.  Termination occurs when the cosine
!      of the angle between fvec and any column of the jacobian is at most
!      gtol in absolute value.  Therefore, gtol measures the orthogonality
!      desired between the function vector and the columns of the jacobian.

!    maxfev is a positive integer input variable.  Termination occurs when the
!      number of calls to fcn is at least maxfev by the end of an iteration.

!    epsfcn is an input variable used in determining a suitable step length
!      for the forward-difference approximation.  This approximation assumes
!      that the relative errors in the functions are of the order of epsfcn.
!      If epsfcn is less than the machine precision, it is assumed that the
!      relative errors in the functions are of the order of the machine
!      precision.

!    diag is an array of length n.  If mode = 1 (see below), diag is
!      internally set.  If mode = 2, diag must contain positive entries that
!      serve as multiplicative scale factors for the variables.

!    mode is an integer input variable.  If mode = 1, the variables will be
!      scaled internally.  If mode = 2, the scaling is specified by the input
!      diag. other values of mode are equivalent to mode = 1.

!    factor is a positive input variable used in determining the initial step
!      bound.  This bound is set to the product of factor and the euclidean
!      norm of diag*x if nonzero, or else to factor itself.  In most cases
!      factor should lie in the interval (.1,100.). 100. is a generally
!      recommended value.

!    nprint is an integer input variable that enables controlled printing of
!      iterates if it is positive.  In this case, fcn is called with iflag = 0
!      at the beginning of the first iteration and every nprint iterations
!      thereafter and immediately prior to return, with x and fvec available
!      for printing.  If nprint is not positive, no special calls
!      of fcn with iflag = 0 are made.

!    info is an integer output variable.  If the user has terminated
!      execution, info is set to the (negative) value of iflag.
!      See description of fcn.  Otherwise, info is set as follows.

!      info = 0  improper input parameters.

!      info = 1  both actual and predicted relative reductions
!                in the sum of squares are at most ftol.

!      info = 2  relative error between two consecutive iterates <= xtol.

!      info = 3  conditions for info = 1 and info = 2 both hold.

!      info = 4  the cosine of the angle between fvec and any column of
!                the Jacobian is at most gtol in absolute value.

!      info = 5  number of calls to fcn has reached or exceeded maxfev.

!      info = 6  ftol is too small. no further reduction in
!                the sum of squares is possible.

!      info = 7  xtol is too small. no further improvement in
!                the approximate solution x is possible.

!      info = 8  gtol is too small. fvec is orthogonal to the
!                columns of the jacobian to machine precision.

!    nfev is an integer output variable set to the number of calls to fcn.

!    fjac is an output m by n array. the upper n by n submatrix
!      of fjac contains an upper triangular matrix r with
!      diagonal elements of nonincreasing magnitude such that

!             t     t           t
!            p *(jac *jac)*p = r *r,

!      where p is a permutation matrix and jac is the final calculated
!      Jacobian.  Column j of p is column ipvt(j) (see below) of the
!      identity matrix. the lower trapezoidal part of fjac contains
!      information generated during the computation of r.

!    ldfjac is a positive integer input variable not less than m
!      which specifies the leading dimension of the array fjac.

!    ipvt is an integer output array of length n.  ipvt defines a permutation
!      matrix p such that jac*p = q*r, where jac is the final calculated
!      jacobian, q is orthogonal (not stored), and r is upper triangular
!      with diagonal elements of nonincreasing magnitude.
!      Column j of p is column ipvt(j) of the identity matrix.

!    qtf is an output array of length n which contains
!      the first n elements of the vector (q transpose)*fvec.

!    wa1, wa2, and wa3 are work arrays of length n.

!    wa4 is a work array of length m.

!  subprograms called

!    user-supplied ...... fcn

!    minpack-supplied ... dpmpar,enorm,fdjac2,lmpar,qrfac

!    fortran-supplied ... dabs,dmax1,dmin1,dsqrt,mod

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  **********
INTEGER   :: i, iflag, iter, j, l
REAL (dp) :: actred, delta, dirder, epsmch, fnorm, fnorm1, gnorm,  &
             par, pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm
REAL (dp) :: diag(n), qtf(n), wa1(n), wa2(n), wa3(n), wa4(m)
REAL (dp), PARAMETER :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp,  &
                        p25 = 0.25_dp, p75 = 0.75_dp, p0001 = 0.0001_dp, &
                        zero = 0.0_dp

!     epsmch is the machine precision.

epsmch = EPSILON(zero)

info = 0
iflag = 0
nfev = 0

!     check the input parameters for errors.

IF (n <= 0 .OR. m < n .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero  &
    .OR. maxfev <= 0 .OR. factor <= zero) GO TO 300
IF (mode /= 2) GO TO 20
DO  j = 1, n
  IF (diag(j) <= zero) GO TO 300
END DO

!     evaluate the function at the starting point and calculate its norm.

20 iflag = 1
CALL fcn(m, n, x, fvec, iflag)
nfev = 1
IF (iflag < 0) GO TO 300
fnorm = enorm(m, fvec)

!     initialize levenberg-marquardt parameter and iteration counter.

par = zero
iter = 1

!     beginning of the outer loop.

!        calculate the jacobian matrix.

30 iflag = 2
CALL fdjac2(fcn, m, n, x, fvec, fjac, iflag, epsfcn)
nfev = nfev + n
IF (iflag < 0) GO TO 300

!        If requested, call fcn to enable printing of iterates.

IF (nprint <= 0) GO TO 40
iflag = 0
IF (MOD(iter-1,nprint) == 0) CALL fcn(m, n, x, fvec, iflag)
IF (iflag < 0) GO TO 300

!        Compute the qr factorization of the jacobian.

40 CALL qrfac(m, n, fjac, .true., ipvt, wa1, wa2)

!        On the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.

IF (iter /= 1) GO TO 80
IF (mode == 2) GO TO 60
DO  j = 1, n
  diag(j) = wa2(j)
  IF (wa2(j) == zero) diag(j) = one
END DO

!        On the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.

60 wa3(1:n) = diag(1:n)*x(1:n)
xnorm = enorm(n, wa3)
delta = factor*xnorm
IF (delta == zero) delta = factor

!        Form (q transpose)*fvec and store the first n components in qtf.

80 wa4(1:m) = fvec(1:m)
DO  j = 1, n
  IF (fjac(j,j) == zero) GO TO 120
  sum = DOT_PRODUCT( fjac(j:m,j), wa4(j:m) )
  temp = -sum/fjac(j,j)
  DO  i = j, m
    wa4(i) = wa4(i) + fjac(i,j)*temp
  END DO
  120 fjac(j,j) = wa1(j)
  qtf(j) = wa4(j)
END DO

!        compute the norm of the scaled gradient.

gnorm = zero
IF (fnorm == zero) GO TO 170
DO  j = 1, n
  l = ipvt(j)
  IF (wa2(l) == zero) CYCLE
  sum = zero
  DO  i = 1, j
    sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  END DO
  gnorm = MAX(gnorm, ABS(sum/wa2(l)))
END DO

!        test for convergence of the gradient norm.

170 IF (gnorm <= gtol) info = 4
IF (info /= 0) GO TO 300

!        rescale if necessary.

IF (mode == 2) GO TO 200
DO  j = 1, n
  diag(j) = MAX(diag(j), wa2(j))
END DO

!        beginning of the inner loop.

!           determine the Levenberg-Marquardt parameter.

200 CALL lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2)

!           store the direction p and x + p. calculate the norm of p.

DO  j = 1, n
  wa1(j) = -wa1(j)
  wa2(j) = x(j) + wa1(j)
  wa3(j) = diag(j)*wa1(j)
END DO
pnorm = enorm(n, wa3)

!           on the first iteration, adjust the initial step bound.

IF (iter == 1) delta = MIN(delta, pnorm)

!           evaluate the function at x + p and calculate its norm.

iflag = 1
CALL fcn(m, n, wa2, wa4, iflag)
nfev = nfev + 1
IF (iflag < 0) GO TO 300
fnorm1 = enorm(m, wa4)

!           compute the scaled actual reduction.

actred = -one
IF (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

!           Compute the scaled predicted reduction and
!           the scaled directional derivative.

DO  j = 1, n
  wa3(j) = zero
  l = ipvt(j)
  temp = wa1(l)
  DO  i = 1, j
    wa3(i) = wa3(i) + fjac(i,j)*temp
  END DO
END DO
temp1 = enorm(n,wa3)/fnorm
temp2 = (SQRT(par)*pnorm)/fnorm
prered = temp1**2 + temp2**2/p5
dirder = -(temp1**2 + temp2**2)

!           compute the ratio of the actual to the predicted reduction.

ratio = zero
IF (prered /= zero) ratio = actred/prered

!           update the step bound.

IF (ratio <= p25) THEN
  IF (actred >= zero) temp = p5
  IF (actred < zero) temp = p5*dirder/(dirder + p5*actred)
  IF (p1*fnorm1 >= fnorm .OR. temp < p1) temp = p1
  delta = temp*MIN(delta,pnorm/p1)
  par = par/temp
ELSE
  IF (par /= zero .AND. ratio < p75) GO TO 260
  delta = pnorm/p5
  par = p5*par
END IF

!           test for successful iteration.

260 IF (ratio < p0001) GO TO 290

!           successful iteration. update x, fvec, and their norms.

DO  j = 1, n
  x(j) = wa2(j)
  wa2(j) = diag(j)*x(j)
END DO
fvec(1:m) = wa4(1:m)
xnorm = enorm(n, wa2)
fnorm = fnorm1
iter = iter + 1

!           tests for convergence.

290 IF (ABS(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one) info = 1
IF (delta <= xtol*xnorm) info = 2
IF (ABS(actred) <= ftol .AND. prered <= ftol  &
    .AND. p5*ratio <= one .AND. info == 2) info = 3
IF (info /= 0) GO TO 300

!           tests for termination and stringent tolerances.

IF (nfev >= maxfev) info = 5
IF (ABS(actred) <= epsmch .AND. prered <= epsmch  &
    .AND. p5*ratio <= one) info = 6
IF (delta <= epsmch*xnorm) info = 7
IF (gnorm <= epsmch) info = 8
IF (info /= 0) GO TO 300

!           end of the inner loop. repeat if iteration unsuccessful.

IF (ratio < p0001) GO TO 200

!        end of the outer loop.

GO TO 30

!     termination, either normal or user imposed.

300 IF (iflag < 0) info = iflag
iflag = 0
IF (nprint > 0) CALL fcn(m, n, x, fvec, iflag)
RETURN

!     last card of subroutine lmdif.

END SUBROUTINE lmdif



SUBROUTINE lmder1(fcn, m, n, x, fvec, fjac, tol, info, ipvt)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:45:54

! N.B. Arguments LDFJAC, WA & LWA have been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(OUT)     :: fvec(:)
REAL (dp), INTENT(IN OUT)  :: fjac(:,:)    ! fjac(ldfjac,n)
REAL (dp), INTENT(IN)      :: tol
INTEGER, INTENT(OUT)       :: info
INTEGER, INTENT(IN OUT)    :: ipvt(:)


! EXTERNAL fcn

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, fjac, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER         :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    REAL (dp), INTENT(OUT)     :: fjac(:,:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn
END INTERFACE

!  **********

!  subroutine lmder1

!  The purpose of lmder1 is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of the
!  levenberg-marquardt algorithm.  This is done by using the more
!  general least-squares solver lmder.  The user must provide a
!  subroutine which calculates the functions and the jacobian.

!  the subroutine statement is

!    subroutine lmder1(fcn, m, n, x, fvec, fjac, tol, info, ipvt)

!  where

!    fcn is the name of the user-supplied subroutine which
!      calculates the functions and the jacobian.  fcn must
!      be declared in an interface statement in the user
!      calling program, and should be written as follows.

!      subroutine fcn(m, n, x, fvec, fjac, iflag)
!      integer   :: m, n, ldfjac, iflag
!      REAL (dp) :: x(:), fvec(:), fjac(:,:)
!      ----------
!      if iflag = 1 calculate the functions at x and
!      return this vector in fvec. do not alter fjac.
!      if iflag = 2 calculate the jacobian at x and
!      return this matrix in fjac. do not alter fvec.
!      ----------
!      return
!      end

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of lmder1.
!      in this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number of functions.

!    n is a positive integer input variable set to the number
!      of variables.  n must not exceed m.

!    x is an array of length n. on input x must contain
!      an initial estimate of the solution vector. on output x
!      contains the final estimate of the solution vector.

!    fvec is an output array of length m which contains
!      the functions evaluated at the output x.

!    fjac is an output m by n array. the upper n by n submatrix
!      of fjac contains an upper triangular matrix r with
!      diagonal elements of nonincreasing magnitude such that

!             t     t           t
!            p *(jac *jac)*p = r *r,

!      where p is a permutation matrix and jac is the final calculated
!      Jacobian.  Column j of p is column ipvt(j) (see below) of the
!      identity matrix.  The lower trapezoidal part of fjac contains
!      information generated during the computation of r.

!    ldfjac is a positive integer input variable not less than m
!      which specifies the leading dimension of the array fjac.

!    tol is a nonnegative input variable. termination occurs
!      when the algorithm estimates either that the relative
!      error in the sum of squares is at most tol or that
!      the relative error between x and the solution is at most tol.

!    info is an integer output variable.  If the user has terminated
!      execution, info is set to the (negative) value of iflag.
!      See description of fcn.  Otherwise, info is set as follows.

!      info = 0  improper input parameters.

!      info = 1  algorithm estimates that the relative error
!                in the sum of squares is at most tol.

!      info = 2  algorithm estimates that the relative error
!                between x and the solution is at most tol.

!      info = 3  conditions for info = 1 and info = 2 both hold.

!      info = 4  fvec is orthogonal to the columns of the
!                jacobian to machine precision.

!      info = 5  number of calls to fcn with iflag = 1 has reached 100*(n+1).

!      info = 6  tol is too small.  No further reduction in
!                the sum of squares is possible.

!      info = 7  tol is too small.  No further improvement in
!                the approximate solution x is possible.

!    ipvt is an integer output array of length n. ipvt
!      defines a permutation matrix p such that jac*p = q*r,
!      where jac is the final calculated jacobian, q is
!      orthogonal (not stored), and r is upper triangular
!      with diagonal elements of nonincreasing magnitude.
!      column j of p is column ipvt(j) of the identity matrix.

!    wa is a work array of length lwa.

!    lwa is a positive integer input variable not less than 5*n+m.

!  subprograms called

!    user-supplied ...... fcn

!    minpack-supplied ... lmder

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  **********
INTEGER   :: maxfev, mode, nfev, njev, nprint
REAL (dp) :: ftol, gtol, xtol
REAL (dp), PARAMETER :: factor = 100._dp, zero = 0.0_dp

info = 0

!     check the input parameters for errors.

IF ( n <= 0 .OR. m < n .OR. tol < zero ) GO TO 10

!     call lmder.

maxfev = 100*(n + 1)
ftol = tol
xtol = tol
gtol = zero
mode = 1
nprint = 0
CALL lmder(fcn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev,  &
           mode, factor, nprint, info, nfev, njev, ipvt)
IF (info == 8) info = 4

10 RETURN

!     last card of subroutine lmder1.

END SUBROUTINE lmder1



SUBROUTINE lmder(fcn, m, n, x, fvec, fjac, ftol, xtol, gtol, maxfev, &
                 mode, factor, nprint, info, nfev, njev, ipvt)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:45:50

! N.B. Arguments LDFJAC, DIAG, QTF, WA1, WA2, WA3 & WA4 have been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(OUT)     :: fvec(m)
REAL (dp), INTENT(OUT)     :: fjac(:,:)    ! fjac(ldfjac,n)
REAL (dp), INTENT(IN)      :: ftol
REAL (dp), INTENT(IN)      :: xtol
REAL (dp), INTENT(IN OUT)  :: gtol
INTEGER, INTENT(IN OUT)    :: maxfev
INTEGER, INTENT(IN)        :: mode
REAL (dp), INTENT(IN)      :: factor
INTEGER, INTENT(IN)        :: nprint
INTEGER, INTENT(OUT)       :: info
INTEGER, INTENT(OUT)       :: nfev
INTEGER, INTENT(OUT)       :: njev
INTEGER, INTENT(OUT)       :: ipvt(:)

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, fjac, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    REAL (dp), INTENT(OUT)     :: fjac(:,:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn
END INTERFACE


!  **********

!  subroutine lmder

!  the purpose of lmder is to minimize the sum of the squares of
!  m nonlinear functions in n variables by a modification of
!  the levenberg-marquardt algorithm. the user must provide a
!  subroutine which calculates the functions and the jacobian.

!  the subroutine statement is

!    subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol,
!                     maxfev,diag,mode,factor,nprint,info,nfev,
!                     njev,ipvt,qtf,wa1,wa2,wa3,wa4)

!  where

!    fcn is the name of the user-supplied subroutine which
!      calculates the functions and the jacobian. fcn must
!      be declared in an external statement in the user
!      calling program, and should be written as follows.

!      subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag)
!      integer m,n,ldfjac,iflag
!      REAL (dp) x(:),fvec(m),fjac(ldfjac,n)
!      ----------
!      if iflag = 1 calculate the functions at x and
!      return this vector in fvec. do not alter fjac.
!      if iflag = 2 calculate the jacobian at x and
!      return this matrix in fjac.  Do not alter fvec.
!      ----------
!      return
!      end

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of lmder.
!      in this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number
!      of functions.

!    n is a positive integer input variable set to the number
!      of variables. n must not exceed m.

!    x is an array of length n. on input x must contain
!      an initial estimate of the solution vector. on output x
!      contains the final estimate of the solution vector.

!    fvec is an output array of length m which contains
!      the functions evaluated at the output x.

!    fjac is an output m by n array. the upper n by n submatrix
!      of fjac contains an upper triangular matrix r with
!      diagonal elements of nonincreasing magnitude such that

!             t     t           t
!            p *(jac *jac)*p = r *r

!      where p is a permutation matrix and jac is the final calculated
!      jacobian.  Column j of p is column ipvt(j) (see below) of the
!      identity matrix.  The lower trapezoidal part of fjac contains
!      information generated during the computation of r.

!    ldfjac is a positive integer input variable not less than m
!      which specifies the leading dimension of the array fjac.

!    ftol is a nonnegative input variable.  Termination occurs when both
!      the actual and predicted relative reductions in the sum of squares
!      are at most ftol.   Therefore, ftol measures the relative error
!      desired in the sum of squares.

!    xtol is a nonnegative input variable. termination
!      occurs when the relative error between two consecutive
!      iterates is at most xtol. therefore, xtol measures the
!      relative error desired in the approximate solution.

!    gtol is a nonnegative input variable.  Termination occurs when the
!      cosine of the angle between fvec and any column of the jacobian is
!      at most gtol in absolute value.  Therefore, gtol measures the
!      orthogonality desired between the function vector and the columns
!      of the jacobian.

!    maxfev is a positive integer input variable.  Termination occurs when
!      the number of calls to fcn with iflag = 1 has reached maxfev.

!    diag is an array of length n.  If mode = 1 (see below), diag is
!      internally set.  If mode = 2, diag must contain positive entries
!      that serve as multiplicative scale factors for the variables.

!    mode is an integer input variable.  if mode = 1, the
!      variables will be scaled internally.  if mode = 2,
!      the scaling is specified by the input diag.  other
!      values of mode are equivalent to mode = 1.

!    factor is a positive input variable used in determining the
!      initial step bound. this bound is set to the product of
!      factor and the euclidean norm of diag*x if nonzero, or else
!      to factor itself. in most cases factor should lie in the
!      interval (.1,100.).100. is a generally recommended value.

!    nprint is an integer input variable that enables controlled printing
!      of iterates if it is positive.  In this case, fcn is called with
!      iflag = 0 at the beginning of the first iteration and every nprint
!      iterations thereafter and immediately prior to return, with x, fvec,
!      and fjac available for printing.  fvec and fjac should not be
!      altered.  If nprint is not positive, no special calls of fcn with
!      iflag = 0 are made.

!    info is an integer output variable.  If the user has terminated
!      execution, info is set to the (negative) value of iflag.
!      See description of fcn.  Otherwise, info is set as follows.

!      info = 0  improper input parameters.

!      info = 1  both actual and predicted relative reductions
!                in the sum of squares are at most ftol.

!      info = 2  relative error between two consecutive iterates
!                is at most xtol.

!      info = 3  conditions for info = 1 and info = 2 both hold.

!      info = 4  the cosine of the angle between fvec and any column of
!                the jacobian is at most gtol in absolute value.

!      info = 5  number of calls to fcn with iflag = 1 has reached maxfev.

!      info = 6  ftol is too small.  No further reduction in
!                the sum of squares is possible.

!      info = 7  xtol is too small.  No further improvement in
!                the approximate solution x is possible.

!      info = 8  gtol is too small.  fvec is orthogonal to the
!                columns of the jacobian to machine precision.

!    nfev is an integer output variable set to the number of
!      calls to fcn with iflag = 1.

!    njev is an integer output variable set to the number of
!      calls to fcn with iflag = 2.

!    ipvt is an integer output array of length n.  ipvt
!      defines a permutation matrix p such that jac*p = q*r,
!      where jac is the final calculated jacobian, q is
!      orthogonal (not stored), and r is upper triangular
!      with diagonal elements of nonincreasing magnitude.
!      column j of p is column ipvt(j) of the identity matrix.

!    qtf is an output array of length n which contains
!      the first n elements of the vector (q transpose)*fvec.

!    wa1, wa2, and wa3 are work arrays of length n.

!    wa4 is a work array of length m.

!  subprograms called

!    user-supplied ...... fcn

!    minpack-supplied ... dpmpar,enorm,lmpar,qrfac

!    fortran-supplied ... ABS,MAX,MIN,SQRT,mod

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  **********
INTEGER   :: i, iflag, iter, j, l
REAL (dp) :: actred, delta, dirder, epsmch, fnorm, fnorm1, gnorm,  &
             par, pnorm, prered, ratio, sum, temp, temp1, temp2, xnorm
REAL (dp) :: diag(n), qtf(n), wa1(n), wa2(n), wa3(n), wa4(m)
REAL (dp), PARAMETER :: one = 1.0_dp, p1 = 0.1_dp, p5 = 0.5_dp,  &
                        p25 = 0.25_dp, p75 = 0.75_dp, p0001 = 0.0001_dp, &
                        zero = 0.0_dp

!     epsmch is the machine precision.

epsmch = EPSILON(zero)

info = 0
iflag = 0
nfev = 0
njev = 0

!     check the input parameters for errors.

IF (n <= 0 .OR. m < n .OR. ftol < zero .OR. xtol < zero .OR. gtol < zero  &
    .OR. maxfev <= 0 .OR. factor <= zero) GO TO 300
IF (mode /= 2) GO TO 20
DO  j = 1, n
  IF (diag(j) <= zero) GO TO 300
END DO

!     evaluate the function at the starting point and calculate its norm.

20 iflag = 1
CALL fcn(m, n, x, fvec, fjac, iflag)
nfev = 1
IF (iflag < 0) GO TO 300
fnorm = enorm(m, fvec)

!     initialize levenberg-marquardt parameter and iteration counter.

par = zero
iter = 1

!     beginning of the outer loop.

!        calculate the jacobian matrix.

30 iflag = 2
CALL fcn(m, n, x, fvec, fjac, iflag)
njev = njev + 1
IF (iflag < 0) GO TO 300

!        if requested, call fcn to enable printing of iterates.

IF (nprint <= 0) GO TO 40
iflag = 0
IF (MOD(iter-1,nprint) == 0) CALL fcn(m, n, x, fvec, fjac, iflag)
IF (iflag < 0) GO TO 300

!        compute the qr factorization of the jacobian.

40 CALL qrfac(m, n, fjac, .true., ipvt, wa1, wa2)

!        on the first iteration and if mode is 1, scale according
!        to the norms of the columns of the initial jacobian.

IF (iter /= 1) GO TO 80
IF (mode == 2) GO TO 60
DO  j = 1, n
  diag(j) = wa2(j)
  IF (wa2(j) == zero) diag(j) = one
END DO

!        on the first iteration, calculate the norm of the scaled x
!        and initialize the step bound delta.

60 wa3(1:n) = diag(1:n)*x(1:n)
xnorm = enorm(n,wa3)
delta = factor*xnorm
IF (delta == zero) delta = factor

!        form (q transpose)*fvec and store the first n components in qtf.

80 wa4(1:m) = fvec(1:m)
DO  j = 1, n
  IF (fjac(j,j) == zero) GO TO 120
  sum = DOT_PRODUCT( fjac(j:m,j), wa4(j:m) )
  temp = -sum/fjac(j,j)
  DO  i = j, m
    wa4(i) = wa4(i) + fjac(i,j)*temp
  END DO
  120 fjac(j,j) = wa1(j)
  qtf(j) = wa4(j)
END DO

!        compute the norm of the scaled gradient.

gnorm = zero
IF (fnorm == zero) GO TO 170
DO  j = 1, n
  l = ipvt(j)
  IF (wa2(l) == zero) CYCLE
  sum = zero
  DO  i = 1, j
    sum = sum + fjac(i,j)*(qtf(i)/fnorm)
  END DO
  gnorm = MAX(gnorm,ABS(sum/wa2(l)))
END DO

!        test for convergence of the gradient norm.

170 IF (gnorm <= gtol) info = 4
IF (info /= 0) GO TO 300

!        rescale if necessary.

IF (mode == 2) GO TO 200
DO  j = 1, n
  diag(j) = MAX(diag(j), wa2(j))
END DO

!        beginning of the inner loop.

!           determine the levenberg-marquardt parameter.

200 CALL lmpar(n, fjac, ipvt, diag, qtf, delta, par, wa1, wa2)

!           store the direction p and x + p. calculate the norm of p.

DO  j = 1, n
  wa1(j) = -wa1(j)
  wa2(j) = x(j) + wa1(j)
  wa3(j) = diag(j)*wa1(j)
END DO
pnorm = enorm(n, wa3)

!           on the first iteration, adjust the initial step bound.

IF (iter == 1) delta = MIN(delta,pnorm)

!           evaluate the function at x + p and calculate its norm.

iflag = 1
CALL fcn(m, n, wa2, wa4, fjac, iflag)
nfev = nfev + 1
IF (iflag < 0) GO TO 300
fnorm1 = enorm(m, wa4)

!           compute the scaled actual reduction.

actred = -one
IF (p1*fnorm1 < fnorm) actred = one - (fnorm1/fnorm)**2

!           compute the scaled predicted reduction and
!           the scaled directional derivative.

DO  j = 1, n
  wa3(j) = zero
  l = ipvt(j)
  temp = wa1(l)
  wa3(1:j) = wa3(1:j) + fjac(1:j,j)*temp
END DO
temp1 = enorm(n,wa3)/fnorm
temp2 = (SQRT(par)*pnorm)/fnorm
prered = temp1**2 + temp2**2/p5
dirder = -(temp1**2 + temp2**2)

!           compute the ratio of the actual to the predicted reduction.

ratio = zero
IF (prered /= zero) ratio = actred/prered

!           update the step bound.

IF (ratio <= p25) THEN
  IF (actred >= zero) temp = p5
  IF (actred < zero) temp = p5*dirder/(dirder + p5*actred)
  IF (p1*fnorm1 >= fnorm .OR. temp < p1) temp = p1
  delta = temp*MIN(delta, pnorm/p1)
  par = par/temp
ELSE
  IF (par /= zero .AND. ratio < p75) GO TO 260
  delta = pnorm/p5
  par = p5*par
END IF

!           test for successful iteration.

260 IF (ratio < p0001) GO TO 290

!           successful iteration. update x, fvec, and their norms.

DO  j = 1, n
  x(j) = wa2(j)
  wa2(j) = diag(j)*x(j)
END DO
fvec(1:m) = wa4(1:m)
xnorm = enorm(n,wa2)
fnorm = fnorm1
iter = iter + 1

!           tests for convergence.

290 IF (ABS(actred) <= ftol .AND. prered <= ftol .AND. p5*ratio <= one) info = 1
IF (delta <= xtol*xnorm) info = 2
IF (ABS(actred) <= ftol .AND. prered <= ftol  &
    .AND. p5*ratio <= one .AND. info == 2) info = 3
IF (info /= 0) GO TO 300

!           tests for termination and stringent tolerances.

IF (nfev >= maxfev) info = 5
IF (ABS(actred) <= epsmch .AND. prered <= epsmch  &
    .AND. p5*ratio <= one) info = 6
IF (delta <= epsmch*xnorm) info = 7
IF (gnorm <= epsmch) info = 8
IF (info /= 0) GO TO 300

!           end of the inner loop. repeat if iteration unsuccessful.

IF (ratio < p0001) GO TO 200

!        end of the outer loop.

GO TO 30

!     termination, either normal or user imposed.

300 IF (iflag < 0) info = iflag
iflag = 0
IF (nprint > 0) CALL fcn(m, n, x, fvec, fjac, iflag)
RETURN

!     last card of subroutine lmder.

END SUBROUTINE lmder



SUBROUTINE lmpar(n, r, ipvt, diag, qtb, delta, par, x, sdiag)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:46:12

! N.B. Arguments LDR, WA1 & WA2 have been removed.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: r(:,:)
INTEGER, INTENT(IN)        :: ipvt(:)
REAL (dp), INTENT(IN)      :: diag(:)
REAL (dp), INTENT(IN)      :: qtb(:)
REAL (dp), INTENT(IN)      :: delta
REAL (dp), INTENT(OUT)     :: par
REAL (dp), INTENT(OUT)     :: x(:)
REAL (dp), INTENT(OUT)     :: sdiag(:)

!  **********

!  subroutine lmpar

!  Given an m by n matrix a, an n by n nonsingular diagonal matrix d,
!  an m-vector b, and a positive number delta, the problem is to determine a
!  value for the parameter par such that if x solves the system

!        a*x = b ,     sqrt(par)*d*x = 0 ,

!  in the least squares sense, and dxnorm is the euclidean
!  norm of d*x, then either par is zero and

!        (dxnorm-delta) <= 0.1*delta ,

!  or par is positive and

!        abs(dxnorm-delta) <= 0.1*delta .

!  This subroutine completes the solution of the problem if it is provided
!  with the necessary information from the r factorization, with column
!  qpivoting, of a.  That is, if a*p = q*r, where p is a permutation matrix,
!  q has orthogonal columns, and r is an upper triangular matrix with diagonal
!  elements of nonincreasing magnitude, then lmpar expects the full upper
!  triangle of r, the permutation matrix p, and the first n components of
!  (q transpose)*b.
!  On output lmpar also provides an upper triangular matrix s such that

!         t   t                   t
!        p *(a *a + par*d*d)*p = s *s .

!  s is employed within lmpar and may be of separate interest.

!  Only a few iterations are generally needed for convergence of the algorithm.
!  If, however, the limit of 10 iterations is reached, then the output par
!  will contain the best value obtained so far.

!  the subroutine statement is

!    subroutine lmpar(n,r,ldr,ipvt,diag,qtb,delta,par,x,sdiag, wa1,wa2)

!  where

!    n is a positive integer input variable set to the order of r.

!    r is an n by n array. on input the full upper triangle
!      must contain the full upper triangle of the matrix r.
!      On output the full upper triangle is unaltered, and the
!      strict lower triangle contains the strict upper triangle
!      (transposed) of the upper triangular matrix s.

!    ldr is a positive integer input variable not less than n
!      which specifies the leading dimension of the array r.

!    ipvt is an integer input array of length n which defines the
!      permutation matrix p such that a*p = q*r. column j of p
!      is column ipvt(j) of the identity matrix.

!    diag is an input array of length n which must contain the
!      diagonal elements of the matrix d.

!    qtb is an input array of length n which must contain the first
!      n elements of the vector (q transpose)*b.

!    delta is a positive input variable which specifies an upper
!      bound on the euclidean norm of d*x.

!    par is a nonnegative variable. on input par contains an
!      initial estimate of the levenberg-marquardt parameter.
!      on output par contains the final estimate.

!    x is an output array of length n which contains the least
!      squares solution of the system a*x = b, sqrt(par)*d*x = 0,
!      for the output par.

!    sdiag is an output array of length n which contains the
!      diagonal elements of the upper triangular matrix s.

!    wa1 and wa2 are work arrays of length n.

!  subprograms called

!    minpack-supplied ... dpmpar,enorm,qrsolv

!    fortran-supplied ... ABS,MAX,MIN,SQRT

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  **********
INTEGER   :: iter, j, jm1, jp1, k, l, nsing
REAL (dp) :: dxnorm, dwarf, fp, gnorm, parc, parl, paru, sum, temp
REAL (dp) :: wa1(n), wa2(n)
REAL (dp), PARAMETER :: p1 = 0.1_dp, p001 = 0.001_dp, zero = 0.0_dp

!     dwarf is the smallest positive magnitude.

dwarf = TINY(zero)

!     compute and store in x the gauss-newton direction. if the
!     jacobian is rank-deficient, obtain a least squares solution.

nsing = n
DO  j = 1, n
  wa1(j) = qtb(j)
  IF (r(j,j) == zero .AND. nsing == n) nsing = j - 1
  IF (nsing < n) wa1(j) = zero
END DO

DO  k = 1, nsing
  j = nsing - k + 1
  wa1(j) = wa1(j)/r(j,j)
  temp = wa1(j)
  jm1 = j - 1
  wa1(1:jm1) = wa1(1:jm1) - r(1:jm1,j)*temp
END DO

DO  j = 1, n
  l = ipvt(j)
  x(l) = wa1(j)
END DO

!     initialize the iteration counter.
!     evaluate the function at the origin, and test
!     for acceptance of the gauss-newton direction.

iter = 0
wa2(1:n) = diag(1:n)*x(1:n)
dxnorm = enorm(n, wa2)
fp = dxnorm - delta
IF (fp <= p1*delta) GO TO 220

!     if the jacobian is not rank deficient, the newton
!     step provides a lower bound, parl, for the zero of
!     the function.  Otherwise set this bound to zero.

parl = zero
IF (nsing < n) GO TO 120
DO  j = 1, n
  l = ipvt(j)
  wa1(j) = diag(l)*(wa2(l)/dxnorm)
END DO
DO  j = 1, n
  sum = DOT_PRODUCT( r(1:j-1,j), wa1(1:j-1) )
  wa1(j) = (wa1(j) - sum)/r(j,j)
END DO
temp = enorm(n,wa1)
parl = ((fp/delta)/temp)/temp

!     calculate an upper bound, paru, for the zero of the function.

120 DO  j = 1, n
  sum = DOT_PRODUCT( r(1:j,j), qtb(1:j) )
  l = ipvt(j)
  wa1(j) = sum/diag(l)
END DO
gnorm = enorm(n,wa1)
paru = gnorm/delta
IF (paru == zero) paru = dwarf/MIN(delta,p1)

!     if the input par lies outside of the interval (parl,paru),
!     set par to the closer endpoint.

par = MAX(par,parl)
par = MIN(par,paru)
IF (par == zero) par = gnorm/dxnorm

!     beginning of an iteration.

150 iter = iter + 1

!        evaluate the function at the current value of par.

IF (par == zero) par = MAX(dwarf, p001*paru)
temp = SQRT(par)
wa1(1:n) = temp*diag(1:n)
CALL qrsolv(n, r, ipvt, wa1, qtb, x, sdiag)
wa2(1:n) = diag(1:n)*x(1:n)
dxnorm = enorm(n, wa2)
temp = fp
fp = dxnorm - delta

!        if the function is small enough, accept the current value
!        of par. also test for the exceptional cases where parl
!        is zero or the number of iterations has reached 10.

IF (ABS(fp) <= p1*delta .OR. parl == zero .AND. fp <= temp  &
    .AND. temp < zero .OR. iter == 10) GO TO 220

!        compute the newton correction.

DO  j = 1, n
  l = ipvt(j)
  wa1(j) = diag(l)*(wa2(l)/dxnorm)
END DO
DO  j = 1, n
  wa1(j) = wa1(j)/sdiag(j)
  temp = wa1(j)
  jp1 = j + 1
  wa1(jp1:n) = wa1(jp1:n) - r(jp1:n,j)*temp
END DO
temp = enorm(n,wa1)
parc = ((fp/delta)/temp)/temp

!        depending on the sign of the function, update parl or paru.

IF (fp > zero) parl = MAX(parl,par)
IF (fp < zero) paru = MIN(paru,par)

!        compute an improved estimate for par.

par = MAX(parl, par+parc)

!        end of an iteration.

GO TO 150

!     termination.

220 IF (iter == 0) par = zero
RETURN

!     last card of subroutine lmpar.

END SUBROUTINE lmpar



SUBROUTINE qrfac(m, n, a, pivot, ipvt, rdiag, acnorm)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:46:17

! N.B. Arguments LDA, LIPVT & WA have been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: a(:,:)
LOGICAL, INTENT(IN)        :: pivot
INTEGER, INTENT(OUT)       :: ipvt(:)
REAL (dp), INTENT(OUT)     :: rdiag(:)
REAL (dp), INTENT(OUT)     :: acnorm(:)

!  **********

!  subroutine qrfac

!  This subroutine uses Householder transformations with column pivoting
!  (optional) to compute a qr factorization of the m by n matrix a.
!  That is, qrfac determines an orthogonal matrix q, a permutation matrix p,
!  and an upper trapezoidal matrix r with diagonal elements of nonincreasing
!  magnitude, such that a*p = q*r.  The householder transformation for
!  column k, k = 1,2,...,min(m,n), is of the form

!                        t
!        i - (1/u(k))*u*u

!  where u has zeros in the first k-1 positions.  The form of this
!  transformation and the method of pivoting first appeared in the
!  corresponding linpack subroutine.

!  the subroutine statement is

!    subroutine qrfac(m, n, a, lda, pivot, ipvt, lipvt, rdiag, acnorm, wa)

! N.B. 3 of these arguments have been omitted in this version.

!  where

!    m is a positive integer input variable set to the number of rows of a.

!    n is a positive integer input variable set to the number of columns of a.

!    a is an m by n array.  On input a contains the matrix for
!      which the qr factorization is to be computed.  On output
!      the strict upper trapezoidal part of a contains the strict
!      upper trapezoidal part of r, and the lower trapezoidal
!      part of a contains a factored form of q (the non-trivial
!      elements of the u vectors described above).

!    lda is a positive integer input variable not less than m
!      which specifies the leading dimension of the array a.

!    pivot is a logical input variable.  If pivot is set true,
!      then column pivoting is enforced.  If pivot is set false,
!      then no column pivoting is done.

!    ipvt is an integer output array of length lipvt.  ipvt
!      defines the permutation matrix p such that a*p = q*r.
!      Column j of p is column ipvt(j) of the identity matrix.
!      If pivot is false, ipvt is not referenced.

!    lipvt is a positive integer input variable.  If pivot is false,
!      then lipvt may be as small as 1.  If pivot is true, then
!      lipvt must be at least n.

!    rdiag is an output array of length n which contains the
!      diagonal elements of r.

!    acnorm is an output array of length n which contains the norms of the
!      corresponding columns of the input matrix a.
!      If this information is not needed, then acnorm can coincide with rdiag.

!    wa is a work array of length n.  If pivot is false, then wa
!      can coincide with rdiag.

!  subprograms called

!    minpack-supplied ... dpmpar,enorm

!    fortran-supplied ... MAX,SQRT,MIN

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  **********
INTEGER   :: i, j, jp1, k, kmax, minmn
REAL (dp) :: ajnorm, epsmch, sum, temp, wa(n)
REAL (dp), PARAMETER :: one = 1.0_dp, p05 = 0.05_dp, zero = 0.0_dp

!     epsmch is the machine precision.

epsmch = EPSILON(zero)

!     compute the initial column norms and initialize several arrays.

DO  j = 1, n
  acnorm(j) = enorm(m,a(1:,j))
  rdiag(j) = acnorm(j)
  wa(j) = rdiag(j)
  IF (pivot) ipvt(j) = j
END DO

!     Reduce a to r with Householder transformations.

minmn = MIN(m,n)
DO  j = 1, minmn
  IF (.NOT.pivot) GO TO 40
  
!        Bring the column of largest norm into the pivot position.
  
  kmax = j
  DO  k = j, n
    IF (rdiag(k) > rdiag(kmax)) kmax = k
  END DO
  IF (kmax == j) GO TO 40
  DO  i = 1, m
    temp = a(i,j)
    a(i,j) = a(i,kmax)
    a(i,kmax) = temp
  END DO
  rdiag(kmax) = rdiag(j)
  wa(kmax) = wa(j)
  k = ipvt(j)
  ipvt(j) = ipvt(kmax)
  ipvt(kmax) = k
  
!     Compute the Householder transformation to reduce the
!     j-th column of a to a multiple of the j-th unit vector.
  
  40 ajnorm = enorm(m-j+1, a(j:,j))
  IF (ajnorm == zero) CYCLE
  IF (a(j,j) < zero) ajnorm = -ajnorm
  a(j:m,j) = a(j:m,j)/ajnorm
  a(j,j) = a(j,j) + one
  
!     Apply the transformation to the remaining columns and update the norms.
  
  jp1 = j + 1
  DO  k = jp1, n
    sum = DOT_PRODUCT( a(j:m,j), a(j:m,k) )
    temp = sum/a(j,j)
    a(j:m,k) = a(j:m,k) - temp*a(j:m,j)
    IF (.NOT.pivot .OR. rdiag(k) == zero) CYCLE
    temp = a(j,k)/rdiag(k)
    rdiag(k) = rdiag(k)*SQRT(MAX(zero, one-temp**2))
    IF (p05*(rdiag(k)/wa(k))**2 > epsmch) CYCLE
    rdiag(k) = enorm(m-j, a(jp1:,k))
    wa(k) = rdiag(k)
  END DO
  rdiag(j) = -ajnorm
END DO

RETURN

!     last card of subroutine qrfac.

END SUBROUTINE qrfac



SUBROUTINE qrsolv(n, r, ipvt, diag, qtb, x, sdiag)
 
! N.B. Arguments LDR & WA have been removed.

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: r(:,:)
INTEGER, INTENT(IN)        :: ipvt(:)
REAL (dp), INTENT(IN)      :: diag(:)
REAL (dp), INTENT(IN)      :: qtb(:)
REAL (dp), INTENT(OUT)     :: x(:)
REAL (dp), INTENT(OUT)     :: sdiag(:)

!  **********

!  subroutine qrsolv

!  Given an m by n matrix a, an n by n diagonal matrix d, and an m-vector b,
!  the problem is to determine an x which solves the system

!        a*x = b ,     d*x = 0 ,

!  in the least squares sense.

!  This subroutine completes the solution of the problem if it is provided
!  with the necessary information from the qr factorization, with column
!  pivoting, of a.  That is, if a*p = q*r, where p is a permutation matrix,
!  q has orthogonal columns, and r is an upper triangular matrix with diagonal
!  elements of nonincreasing magnitude, then qrsolv expects the full upper
!  triangle of r, the permutation matrix p, and the first n components of
!  (q transpose)*b.  The system a*x = b, d*x = 0, is then equivalent to

!               t       t
!        r*z = q *b ,  p *d*p*z = 0 ,

!  where x = p*z. if this system does not have full rank,
!  then a least squares solution is obtained.  On output qrsolv
!  also provides an upper triangular matrix s such that

!         t   t               t
!        p *(a *a + d*d)*p = s *s .

!  s is computed within qrsolv and may be of separate interest.

!  the subroutine statement is

!    subroutine qrsolv(n, r, ldr, ipvt, diag, qtb, x, sdiag, wa)

! N.B. Arguments LDR and WA have been removed in this version.

!  where

!    n is a positive integer input variable set to the order of r.

!    r is an n by n array.  On input the full upper triangle must contain
!      the full upper triangle of the matrix r.
!      On output the full upper triangle is unaltered, and the strict lower
!      triangle contains the strict upper triangle (transposed) of the
!      upper triangular matrix s.

!    ldr is a positive integer input variable not less than n
!      which specifies the leading dimension of the array r.

!    ipvt is an integer input array of length n which defines the
!      permutation matrix p such that a*p = q*r.  Column j of p
!      is column ipvt(j) of the identity matrix.

!    diag is an input array of length n which must contain the
!      diagonal elements of the matrix d.

!    qtb is an input array of length n which must contain the first
!      n elements of the vector (q transpose)*b.

!    x is an output array of length n which contains the least
!      squares solution of the system a*x = b, d*x = 0.

!    sdiag is an output array of length n which contains the
!      diagonal elements of the upper triangular matrix s.

!    wa is a work array of length n.

!  subprograms called

!    fortran-supplied ... ABS,SQRT

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  **********
INTEGER   :: i, j, k, kp1, l, nsing
REAL (dp) :: COS, cotan, qtbpj, SIN, sum, TAN, temp, wa(n)
REAL (dp), PARAMETER :: p5 = 0.5_dp, p25 = 0.25_dp, zero = 0.0_dp

!     Copy r and (q transpose)*b to preserve input and initialize s.
!     In particular, save the diagonal elements of r in x.

DO  j = 1, n
  r(j:n,j) = r(j,j:n)
  x(j) = r(j,j)
  wa(j) = qtb(j)
END DO

!     Eliminate the diagonal matrix d using a givens rotation.

DO  j = 1, n
  
!        Prepare the row of d to be eliminated, locating the
!        diagonal element using p from the qr factorization.
  
  l = ipvt(j)
  IF (diag(l) == zero) CYCLE
  sdiag(j:n) = zero
  sdiag(j) = diag(l)
  
!     The transformations to eliminate the row of d modify only a single
!     element of (q transpose)*b beyond the first n, which is initially zero.
  
  qtbpj = zero
  DO  k = j, n
    
!        Determine a givens rotation which eliminates the
!        appropriate element in the current row of d.
    
    IF (sdiag(k) == zero) CYCLE
    IF (ABS(r(k,k)) < ABS(sdiag(k))) THEN
      cotan = r(k,k)/sdiag(k)
      SIN = p5/SQRT(p25 + p25*cotan**2)
      COS = SIN*cotan
    ELSE
      TAN = sdiag(k)/r(k,k)
      COS = p5/SQRT(p25 + p25*TAN**2)
      SIN = COS*TAN
    END IF
    
!        Compute the modified diagonal element of r and
!        the modified element of ((q transpose)*b,0).
    
    r(k,k) = COS*r(k,k) + SIN*sdiag(k)
    temp = COS*wa(k) + SIN*qtbpj
    qtbpj = -SIN*wa(k) + COS*qtbpj
    wa(k) = temp
    
!        Accumulate the tranformation in the row of s.
    
    kp1 = k + 1
    DO  i = kp1, n
      temp = COS*r(i,k) + SIN*sdiag(i)
      sdiag(i) = -SIN*r(i,k) + COS*sdiag(i)
      r(i,k) = temp
    END DO
  END DO
  
!     Store the diagonal element of s and restore
!     the corresponding diagonal element of r.
  
  sdiag(j) = r(j,j)
  r(j,j) = x(j)
END DO

!     Solve the triangular system for z.  If the system is singular,
!     then obtain a least squares solution.

nsing = n
DO  j = 1, n
  IF (sdiag(j) == zero .AND. nsing == n) nsing = j - 1
  IF (nsing < n) wa(j) = zero
END DO

DO  k = 1, nsing
  j = nsing - k + 1
  sum = DOT_PRODUCT( r(j+1:nsing,j), wa(j+1:nsing) )
  wa(j) = (wa(j) - sum)/sdiag(j)
END DO

!     Permute the components of z back to components of x.

DO  j = 1, n
  l = ipvt(j)
  x(l) = wa(j)
END DO
RETURN

!     last card of subroutine qrsolv.

END SUBROUTINE qrsolv



FUNCTION enorm(n,x) RESULT(fn_val)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:45:34

INTEGER, INTENT(IN)    :: n
REAL (dp), INTENT(IN)  :: x(:)
REAL (dp)              :: fn_val

!  **********

!  function enorm

!  given an n-vector x, this function calculates the euclidean norm of x.

!  the euclidean norm is computed by accumulating the sum of squares in
!  three different sums.  The sums of squares for the small and large
!  components are scaled so that no overflows occur.  Non-destructive
!  underflows are permitted.  Underflows and overflows do not occur in the
!  computation of the unscaled sum of squares for the intermediate
!  components.  The definitions of small, intermediate and large components
!  depend on two constants, rdwarf and rgiant.  The main restrictions on
!  these constants are that rdwarf**2 not underflow and rgiant**2 not
!  overflow.  The constants given here are suitable for every known computer.

!  the function statement is

!    REAL (dp) function enorm(n,x)

!  where

!    n is a positive integer input variable.

!    x is an input array of length n.

!  subprograms called

!    fortran-supplied ... ABS,SQRT

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  **********
INTEGER   :: i
REAL (dp) :: agiant, floatn, s1, s2, s3, xabs, x1max, x3max
REAL (dp), PARAMETER :: one = 1.0_dp, zero = 0.0_dp, rdwarf = 3.834E-20_dp,  &
                        rgiant = 1.304E+19_dp

s1 = zero
s2 = zero
s3 = zero
x1max = zero
x3max = zero
floatn = n
agiant = rgiant/floatn
DO  i = 1, n
  xabs = ABS(x(i))
  IF (xabs > rdwarf .AND. xabs < agiant) GO TO 70
  IF (xabs <= rdwarf) GO TO 30
  
!              sum for large components.
  
  IF (xabs <= x1max) GO TO 10
  s1 = one + s1*(x1max/xabs)**2
  x1max = xabs
  GO TO 20

  10 s1 = s1 + (xabs/x1max)**2

  20 GO TO 60
  
!              sum for small components.
  
  30 IF (xabs <= x3max) GO TO 40
  s3 = one + s3*(x3max/xabs)**2
  x3max = xabs
  GO TO 60

  40 IF (xabs /= zero) s3 = s3 + (xabs/x3max)**2

  60 CYCLE
  
!           sum for intermediate components.
  
  70 s2 = s2 + xabs**2
END DO

!     calculation of norm.

IF (s1 == zero) GO TO 100
fn_val = x1max*SQRT(s1 + (s2/x1max)/x1max)
GO TO 120

100 IF (s2 == zero) GO TO 110
IF (s2 >= x3max) fn_val = SQRT(s2*(one + (x3max/s2)*(x3max*s3)))
IF (s2 < x3max) fn_val = SQRT(x3max*((s2/x3max) + (x3max*s3)))
GO TO 120

110 fn_val = x3max*SQRT(s3)

120 RETURN

!     last card of function enorm.

END FUNCTION enorm



SUBROUTINE fdjac2(fcn, m, n, x, fvec, fjac, iflag, epsfcn)
 
! Code converted using TO_F90 by Alan Miller
! Date: 1999-12-09  Time: 12:45:44

! N.B. Arguments LDFJAC & WA have been removed.

INTEGER, INTENT(IN)        :: m
INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(n)
REAL (dp), INTENT(IN)      :: fvec(m)
REAL (dp), INTENT(OUT)     :: fjac(:,:)    ! fjac(ldfjac,n)
INTEGER, INTENT(IN OUT)    :: iflag
REAL (dp), INTENT(IN)      :: epsfcn

INTERFACE
  SUBROUTINE fcn(m, n, x, fvec, iflag)
    IMPLICIT NONE
    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)
    INTEGER, INTENT(IN)        :: m, n
    REAL (dp), INTENT(IN)      :: x(:)
    REAL (dp), INTENT(IN OUT)  :: fvec(:)
    INTEGER, INTENT(IN OUT)    :: iflag
  END SUBROUTINE fcn
END INTERFACE

!  **********

!  subroutine fdjac2

!  this subroutine computes a forward-difference approximation
!  to the m by n jacobian matrix associated with a specified
!  problem of m functions in n variables.

!  the subroutine statement is

!    subroutine fdjac2(fcn,m,n,x,fvec,fjac,ldfjac,iflag,epsfcn,wa)

!  where

!    fcn is the name of the user-supplied subroutine which calculates the
!      functions.  fcn must be declared in an external statement in the user
!      calling program, and should be written as follows.

!      subroutine fcn(m,n,x,fvec,iflag)
!      integer m,n,iflag
!      REAL (dp) x(n),fvec(m)
!      ----------
!      calculate the functions at x and
!      return this vector in fvec.
!      ----------
!      return
!      end

!      the value of iflag should not be changed by fcn unless
!      the user wants to terminate execution of fdjac2.
!      in this case set iflag to a negative integer.

!    m is a positive integer input variable set to the number of functions.

!    n is a positive integer input variable set to the number of variables.
!      n must not exceed m.

!    x is an input array of length n.

!    fvec is an input array of length m which must contain the
!      functions evaluated at x.

!    fjac is an output m by n array which contains the
!      approximation to the jacobian matrix evaluated at x.

!    ldfjac is a positive integer input variable not less than m
!      which specifies the leading dimension of the array fjac.

!    iflag is an integer variable which can be used to terminate
!      the execution of fdjac2.  see description of fcn.

!    epsfcn is an input variable used in determining a suitable step length
!      for the forward-difference approximation.  This approximation assumes
!      that the relative errors in the functions are of the order of epsfcn.
!      If epsfcn is less than the machine precision, it is assumed that the
!      relative errors in the functions are of the order of the machine
!      precision.

!    wa is a work array of length m.

!  subprograms called

!    user-supplied ...... fcn

!    minpack-supplied ... dpmpar

!    fortran-supplied ... ABS,MAX,SQRT

!  argonne national laboratory. minpack project. march 1980.
!  burton s. garbow, kenneth e. hillstrom, jorge j. more

!  **********
INTEGER   :: j
REAL (dp) :: eps, epsmch, h, temp, wa(m)
REAL (dp), PARAMETER :: zero = 0.0_dp

!     epsmch is the machine precision.

epsmch = EPSILON(zero)

eps = SQRT(MAX(epsfcn, epsmch))
DO  j = 1, n
  temp = x(j)
  h = eps*ABS(temp)
  IF (h == zero) h = eps
  x(j) = temp + h
  CALL fcn(m, n, x, wa, iflag)
  IF (iflag < 0) EXIT
  x(j) = temp
  fjac(1:m,j) = (wa(1:m) - fvec(1:m))/h
END DO

RETURN

!     last card of subroutine fdjac2.

END SUBROUTINE fdjac2

END MODULE Levenberg_Marquardt


module bspline
    
!*****************************************************************************************
!> author: Jacob Williams
!  license: BSD
!
!# Description
!
!  Multidimensional (1D-6D) B-spline interpolation of data on a regular grid.
!  Basic pure subroutine interface.
!
!# Notes
!
!  This module is based on the B-spline and spline routines from [1].
!  The original Fortran 77 routines were converted to free-form source.
!  Some of them are relatively unchanged from the originals, but some have
!  been extensively refactored. In addition, new routines for
!  1d, 4d, 5d, and 6d interpolation were also created (these are simply
!  extensions of the same algorithm into higher dimensions).
!
!# See also
!  * An object-oriented interface can be found in [[bspline_oo_module]].
!
!# References
!
!  1. DBSPLIN and DTENSBS from the
!     [NIST Core Math Library](http://www.nist.gov/itl/math/mcsd-software.cfm).
!     Original code is public domain.
!  2. Carl de Boor, "A Practical Guide to Splines",
!     Springer-Verlag, New York, 1978.
!  3. Carl de Boor, [Efficient Computer Manipulation of Tensor
!     Products](http://dl.acm.org/citation.cfm?id=355831),
!     ACM Transactions on Mathematical Software,
!     Vol. 5 (1979), p. 173-182.
!  4. D.E. Amos, "Computation with Splines and B-Splines",
!     SAND78-1968, Sandia Laboratories, March, 1979.
!  5. Carl de Boor,
!     [Package for calculating with B-splines](http://epubs.siam.org/doi/abs/10.1137/0714026),
!     SIAM Journal on Numerical Analysis 14, 3 (June 1977), p. 441-472.


    use emission, only : dp
    use,intrinsic :: iso_fortran_env, only: error_unit

    implicit none

    private

    integer,parameter :: wp = dp  !! Real precision

    public :: db3ink, db3val

    public :: get_status_message

    contains

!*****************************************************************************************

!*****************************************************************************************
!> Determines the parameters of a function that interpolates
!  the three-dimensional gridded data
!  $$ [x(i),y(j),z(k),\mathrm{fcn}(i,j,k)] ~\mathrm{for}~
!     i=1,..,n_x ~\mathrm{and}~ j=1,..,n_y, ~\mathrm{and}~ k=1,..,n_z $$
!  The interpolating function and
!  its derivatives may subsequently be evaluated by the function
!  [[db3val]].
!
!  The interpolating function is a piecewise polynomial function
!  represented as a tensor product of one-dimensional b-splines. the
!  form of this function is
!  $$ s(x,y,z) = \sum_{i=1}^{n_x} \sum_{j=1}^{n_y} \sum_{k=1}^{n_z}
!                a_{ijk} u_i(x) v_j(y) w_k(z) $$
!
!  where the functions \(u_i\), \(v_j\), and \(w_k\) are one-dimensional b-
!  spline basis functions. the coefficients \(a_{ijk}\) are chosen so that:
!
!  $$ s(x(i),y(j),z(k)) = \mathrm{fcn}(i,j,k)
!     ~\mathrm{for}~ i=1,..,n_x , j=1,..,n_y , k=1,..,n_z $$
!
!  Note that for fixed values of y and z s(x,y,z) is a piecewise
!  polynomial function of x alone, for fixed values of x and z s(x,y,z)
!  is a piecewise polynomial function of y alone, and for fixed
!  values of x and y s(x,y,z) is a function of z alone. in one
!  dimension a piecewise polynomial may be created by partitioning a
!  given interval into subintervals and defining a distinct polynomial
!  piece on each one. the points where adjacent subintervals meet are
!  called knots. each of the functions \(u_i\), \(v_j\), and \(w_k\) above is a
!  piecewise polynomial.
!
!  Users of db3ink choose the order (degree+1) of the polynomial
!  pieces used to define the piecewise polynomial in each of the x, y,
!  and z directions (kx, ky, and kz). users also may define their own
!  knot sequence in x, y, and z separately (tx, ty, and tz). if iflag=
!  0, however, db3ink will choose sequences of knots that result in a
!  piecewise polynomial interpolant with kx-2 continuous partial
!  derivatives in x, ky-2 continuous partial derivatives in y, and kz-
!  2 continuous partial derivatives in z. (kx knots are taken near
!  each endpoint in x, not-a-knot end conditions are used, and the
!  remaining knots are placed at data points if kx is even or at
!  midpoints between data points if kx is odd. the y and z directions
!  are treated similarly.)
!
!  After a call to db3ink, all information necessary to define the
!  interpolating function are contained in the parameters nx, ny, nz,
!  kx, ky, kz, tx, ty, tz, and bcoef. these quantities should not be
!  altered until after the last call of the evaluation routine [[db3val]].
!
!# History
!
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db3ink(x,nx,y,ny,z,nz,fcn,kx,ky,kz,iknot,tx,ty,tz,bcoef,iflag)

    implicit none

    integer,intent(in)                       :: nx    !! number of x abcissae (>= 3)
    integer,intent(in)                       :: ny    !! number of y abcissae (>= 3)
    integer,intent(in)                       :: nz    !! number of z abcissae (>= 3)
    integer,intent(in)                       :: kx    !! the order of spline pieces in x (>= 2, < nx). (order = polynomial degree + 1)
    integer,intent(in)                       :: ky    !! the order of spline pieces in y (>= 2, < ny). (order = polynomial degree + 1)
    integer,intent(in)                       :: kz    !! the order of spline pieces in z (>= 2, < nz). (order = polynomial degree + 1)
    real(wp),dimension(:),intent(in)         :: x     !! `nx` array of x abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)         :: y     !! `ny` array of y abcissae. must be strictly increasing.
    real(wp),dimension(:),intent(in)         :: z     !! `nz` array of z abcissae. must be strictly increasing.
    real(wp),dimension(:,:,:),intent(in)     :: fcn   !! `(nx,ny,nz)` matrix of function values to interpolate. fcn(i,j,k) should
                                                      !!   contain the function value at the point (x(i),y(j),z(k))
    integer,intent(in)                       :: iknot !! 0 = knot sequence chosen by [[db1ink]].
                                                      !! 1 = knot sequence chosen by user.
    real(wp),dimension(:),intent(inout)      :: tx    !! The `nx+kx` knots in the `x` direction for the spline interpolant.
                                                      !!   If `iknot=0` these are chosen by [[db3ink]].
                                                      !!   If `iknot=1` these are specified by the user.
                                                      !!    Must be non-decreasing.
    real(wp),dimension(:),intent(inout)      :: ty    !! The `ny+ky` knots in the `y` direction for the spline interpolant.
                                                      !!    If `iknot=0` these are chosen by [[db3ink]].
                                                      !!    If `iknot=1` these are specified by the user.
                                                      !!    Must be non-decreasing.
    real(wp),dimension(:),intent(inout)      :: tz    !! The `nz+kz` knots in the `z` direction for the spline interpolant.
                                                      !!    If `iknot=0` these are chosen by [[db3ink]].
                                                      !!    If `iknot=1` these are specified by the user.
                                                      !!    Must be non-decreasing.
    real(wp),dimension(:,:,:),intent(out)    :: bcoef !! '(nx,ny,nz)' matrix of coefficients of the b-spline interpolant.
    integer,intent(out)                      :: iflag !!  0 = successful execution.
                                                      !!  2 = iknot out of range.
                                                      !!  3 = nx out of range.
                                                      !!  4 = kx out of range.
                                                      !!  5 = x not strictly increasing.
                                                      !!  6 = tx not non-decreasing.
                                                      !!  7 = ny out of range.
                                                      !!  8 = ky out of range.
                                                      !!  9 = y not strictly increasing.
                                                      !! 10 = ty not non-decreasing.
                                                      !! 11 = nz out of range.
                                                      !! 12 = kz out of range.
                                                      !! 13 = z not strictly increasing.
                                                      !! 14 = ty not non-decreasing.
                                                      !! 700 = size(x) /= size(fcn,1).
                                                      !! 701 = size(y) /= size(fcn,2).
                                                      !! 702 = size(z) /= size(fcn,3).
                                                      !! 706 = size(x) /= nx.
                                                      !! 707 = size(y) /= ny.
                                                      !! 708 = size(z) /= nz.
                                                      !! 712 = size(tx) /= nx+kx.
                                                      !! 713 = size(ty) /= ny+ky.
                                                      !! 714 = size(tz) /= nz+kz.
                                                      !! 800 = size(x) /= size(bcoef,1).
                                                      !! 801 = size(y) /= size(bcoef,2).
                                                      !! 802 = size(z) /= size(bcoef,3).

    real(wp),dimension(nx*ny*nz) :: temp
    real(wp),dimension(max(2*kx*(nx+1),2*ky*(ny+1),2*kz*(nz+1))) :: work
    logical :: status_ok

    ! check validity of input

    call check_inputs('db3ink',&
                        iknot,&
                        iflag,&
                        nx=nx,ny=ny,nz=nz,&
                        kx=kx,ky=ky,kz=kz,&
                        x=x,y=y,z=z,&
                        tx=tx,ty=ty,tz=tz,&
                        f3=fcn,&
                        bcoef3=bcoef,&
                        status_ok=status_ok)

    if (status_ok) then

        ! choose knots

        if (iknot == 0) then
            call dbknot(x,nx,kx,tx)
            call dbknot(y,ny,ky,ty)
            call dbknot(z,nz,kz,tz)
        end if

        ! copy fcn to work in packed for dbtpcf
        temp(1:nx*ny*nz) = reshape( fcn, [nx*ny*nz] )

        ! construct b-spline coefficients

                      call dbtpcf(x,nx,temp, nx,ny*nz,tx,kx,bcoef,work,iflag)
        if (iflag==0) call dbtpcf(y,ny,bcoef,ny,nx*nz,ty,ky,temp, work,iflag)
        if (iflag==0) call dbtpcf(z,nz,temp, nz,nx*ny,tz,kz,bcoef,work,iflag)

    end if

    end subroutine db3ink
!*****************************************************************************************

!*****************************************************************************************
!> Evaluates the tensor product piecewise polynomial
!  interpolant constructed by the routine [[db3ink]] or one of its
!  derivatives at the point (xval,yval,zval).
!
!  To evaluate the
!  interpolant itself, set idx=idy=idz=0, to evaluate the first
!  partial with respect to x, set idx=1,idy=idz=0, and so on.
!
!  db3val returns 0.0 if (xval,yval,zval) is out of range. that is,
!```fortran
! xval<tx(1) .or. xval>tx(nx+kx) .or.
! yval<ty(1) .or. yval>ty(ny+ky) .or.
! zval<tz(1) .or. zval>tz(nz+kz)
!```
!  if the knots tx, ty, and tz were chosen by [[db3ink]], then this is
!  equivalent to
!```fortran
! xval<x(1) .or. xval>x(nx)+epsx .or.
! yval<y(1) .or. yval>y(ny)+epsy .or.
! zval<z(1) .or. zval>z(nz)+epsz
!```
!  where
!```fortran
! epsx = 0.1*(x(nx)-x(nx-1))
! epsy = 0.1*(y(ny)-y(ny-1))
! epsz = 0.1*(z(nz)-z(nz-1))
!```
!
!  The input quantities tx, ty, tz, nx, ny, nz, kx, ky, kz, and bcoef
!  should remain unchanged since the last call of [[db3ink]].
!
!# History
!
!  * Boisvert, Ronald, NBS : 25 may 1982 : Author of original routine.
!  * JEC : 000330 modified array declarations.
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine db3val(xval,yval,zval,idx,idy,idz,&
                                     tx,ty,tz,&
                                     nx,ny,nz,kx,ky,kz,bcoef,f,iflag,&
                                     inbvx,inbvy,inbvz,iloy,iloz)

    implicit none

    integer,intent(in)                      :: idx      !! x derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: idy      !! y derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: idz      !! z derivative of piecewise polynomial to evaluate.
    integer,intent(in)                      :: nx       !! the number of interpolation points in x. (same as in last call to [[db3ink]])
    integer,intent(in)                      :: ny       !! the number of interpolation points in y. (same as in last call to [[db3ink]])
    integer,intent(in)                      :: nz       !! the number of interpolation points in z. (same as in last call to [[db3ink]])
    integer,intent(in)                      :: kx       !! order of polynomial pieces in x. (same as in last call to [[db3ink]])
    integer,intent(in)                      :: ky       !! order of polynomial pieces in y. (same as in last call to [[db3ink]])
    integer,intent(in)                      :: kz       !! order of polynomial pieces in z. (same as in last call to [[db3ink]])
    real(wp),intent(in)                     :: xval     !! x coordinate of evaluation point.
    real(wp),intent(in)                     :: yval     !! y coordinate of evaluation point.
    real(wp),intent(in)                     :: zval     !! z coordinate of evaluation point.
    real(wp),dimension(nx+kx),intent(in)    :: tx       !! sequence of knots defining the piecewise polynomial in the x direction. (same as in last call to [[db3ink]])
    real(wp),dimension(ny+ky),intent(in)    :: ty       !! sequence of knots defining the piecewise polynomial in the y direction. (same as in last call to [[db3ink]])
    real(wp),dimension(nz+kz),intent(in)    :: tz       !! sequence of knots defining the piecewise polynomial in the z direction. (same as in last call to [[db3ink]])
    real(wp),dimension(nx,ny,nz),intent(in) :: bcoef    !! the b-spline coefficients computed by [[db3ink]].
    real(wp),intent(out)                    :: f        !! interpolated value
    integer,intent(out)                     :: iflag    !! status flag: 0 : no errors, /=0 : error
    integer,intent(inout)                   :: inbvx    !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.
    integer,intent(inout)                   :: inbvy    !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.
    integer,intent(inout)                   :: inbvz    !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.
    integer,intent(inout)                   :: iloy     !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.
    integer,intent(inout)                   :: iloz     !! initialization parameter which must be set to 1 the first time this routine is called, and must not be changed by the user.

    real(wp),dimension(ky,kz)              :: temp1
    real(wp),dimension(kz)                 :: temp2
    real(wp),dimension(3*max(kx,ky,kz))    :: work

    integer :: lefty, leftz, mflag,&
                kcoly, kcolz, j, k

    f = 0.0_wp

    if (xval<tx(1) .or. xval>tx(nx+kx)) then
        !write(error_unit,'(A)') 'db3val - x value out of bounds'
        iflag = 601
        return
    end if
    if (yval<ty(1) .or. yval>ty(ny+ky)) then
        !write(error_unit,'(A)') 'db3val - y value out of bounds'
        iflag = 602
        return
    end if
    if (zval<tz(1) .or. zval>tz(nz+kz)) then
        !write(error_unit,'(A)') 'db3val - z value out of bounds'
        iflag = 603
        return
    end if

    iflag = -1
    call dintrv(ty,ny+ky,yval,iloy,lefty,mflag); if (mflag /= 0) return
    call dintrv(tz,nz+kz,zval,iloz,leftz,mflag); if (mflag /= 0) return

    iflag = 0

    kcolz = leftz - kz
    do k=1,kz
        kcolz = kcolz + 1
        kcoly = lefty - ky
        do j=1,ky
            kcoly = kcoly + 1
            call dbvalu(tx,bcoef(:,kcoly,kcolz),nx,kx,idx,xval,inbvx,work,iflag,temp1(j,k))
            if (iflag/=0) return
        end do
    end do

    kcoly = lefty - ky + 1
    do k=1,kz
        call dbvalu(ty(kcoly:),temp1(:,k),ky,ky,idy,yval,inbvy,work,iflag,temp2(k))
        if (iflag/=0) return
    end do

    kcolz = leftz - kz + 1
    call dbvalu(tz(kcolz:),temp2,kz,kz,idz,zval,inbvz,work,iflag,f)

    end subroutine db3val
!*****************************************************************************************

!*****************************************************************************************
!> Check the validity of the inputs to the "ink" routines.
!  Prints warning message if there is an error,
!  and also sets iflag and status_ok.
!
!  Supports up to 6D: x,y,z,q,r,s
!
!# Notes
!
!  The code is new, but the logic is based on the original
!  logic in the CMLIB routines db2ink and db3ink.
!
!# History
!
!  * Jacob Williams, 2/24/2015 : Created this routine.

    pure subroutine check_inputs(routine,&
                            iknot,&
                            iflag,&
                            nx,ny,nz,nq,nr,ns,&
                            kx,ky,kz,kq,kr,ks,&
                            x,y,z,q,r,s,&
                            tx,ty,tz,tq,tr,ts,&
                            f1,f2,f3,f4,f5,f6,&
                            bcoef1,bcoef2,bcoef3,bcoef4,bcoef5,bcoef6,&
                            status_ok)

    implicit none

    character(len=*),intent(in)                         :: routine
    integer,intent(in)                                  :: iknot !! = 0 if the INK routine is computing the knots.
    integer,intent(out)                                 :: iflag
    integer,intent(in),optional                         :: nx,ny,nz,nq,nr,ns
    integer,intent(in),optional                         :: kx,ky,kz,kq,kr,ks
    real(wp),dimension(:),intent(in),optional           :: x,y,z,q,r,s
    real(wp),dimension(:),intent(in),optional           :: tx,ty,tz,tq,tr,ts
    real(wp),dimension(:),intent(in),optional           :: f1,bcoef1
    real(wp),dimension(:,:),intent(in),optional         :: f2,bcoef2
    real(wp),dimension(:,:,:),intent(in),optional       :: f3,bcoef3
    real(wp),dimension(:,:,:,:),intent(in),optional     :: f4,bcoef4
    real(wp),dimension(:,:,:,:,:),intent(in),optional   :: f5,bcoef5
    real(wp),dimension(:,:,:,:,:,:),intent(in),optional :: f6,bcoef6
    logical,intent(out)                                 :: status_ok

    logical :: error

    status_ok = .false.

    if ((iknot < 0) .or. (iknot > 1)) then

        !write(error_unit,'(A,1X,I5)') &
        !    trim(routine)//' - iknot is out of range: ',iflag
        iflag = 2

    else

        call check('x',nx,kx,x,tx,[3,  4, 5, 6,706,712],iflag,error); if (error) return
        call check('y',ny,ky,y,ty,[7,  8, 9,10,707,713],iflag,error); if (error) return
        call check('z',nz,kz,z,tz,[11,12,13,14,708,714],iflag,error); if (error) return
        call check('q',nq,kq,q,tq,[15,16,17,18,709,715],iflag,error); if (error) return
        call check('r',nr,kr,r,tr,[19,20,21,22,710,716],iflag,error); if (error) return
        call check('s',ns,ks,s,ts,[23,24,25,26,711,717],iflag,error); if (error) return

        if (present(x) .and. present(f1) .and. present(bcoef1)) then
            if (size(x)/=size(f1,1)) then; iflag = 700; return; end if
            if (size(x)/=size(bcoef1,1)) then; iflag = 800; return; end if
        end if
        if (present(x) .and. present(y) .and. present(f2) .and. present(bcoef2)) then
            if (size(x)/=size(f2,1)) then; iflag = 700; return; end if
            if (size(y)/=size(f2,2)) then; iflag = 701; return; end if
            if (size(x)/=size(bcoef2,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef2,2)) then; iflag = 801; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(f3) .and. &
            present(bcoef3)) then
            if (size(x)/=size(f3,1)) then; iflag = 700; return; end if
            if (size(y)/=size(f3,2)) then; iflag = 701; return; end if
            if (size(z)/=size(f3,2)) then; iflag = 702; return; end if
            if (size(x)/=size(bcoef3,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef3,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef3,3)) then; iflag = 802; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(f4) .and. present(bcoef4)) then
            if (size(x)/=size(f4,1)) then; iflag = 700; return; end if
            if (size(y)/=size(f4,2)) then; iflag = 701; return; end if
            if (size(z)/=size(f4,3)) then; iflag = 702; return; end if
            if (size(q)/=size(f4,4)) then; iflag = 703; return; end if
            if (size(x)/=size(bcoef4,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef4,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef4,3)) then; iflag = 802; return; end if
            if (size(q)/=size(bcoef4,4)) then; iflag = 803; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(f5) .and. present(bcoef5)) then
            if (size(x)/=size(f5,1)) then; iflag = 700; return; end if
            if (size(y)/=size(f5,2)) then; iflag = 701; return; end if
            if (size(z)/=size(f5,3)) then; iflag = 702; return; end if
            if (size(q)/=size(f5,4)) then; iflag = 703; return; end if
            if (size(r)/=size(f5,5)) then; iflag = 704; return; end if
            if (size(x)/=size(bcoef5,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef5,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef5,3)) then; iflag = 802; return; end if
            if (size(q)/=size(bcoef5,4)) then; iflag = 803; return; end if
            if (size(r)/=size(bcoef5,5)) then; iflag = 804; return; end if
        end if
        if (present(x) .and. present(y) .and. present(z) .and. present(q) .and. &
            present(r) .and. present(s) .and. present(f6) .and. present(bcoef6)) then
            if (size(x)/=size(f6,1)) then; iflag = 700; return; end if
            if (size(y)/=size(f6,2)) then; iflag = 701; return; end if
            if (size(z)/=size(f6,3)) then; iflag = 702; return; end if
            if (size(q)/=size(f6,4)) then; iflag = 703; return; end if
            if (size(r)/=size(f6,5)) then; iflag = 704; return; end if
            if (size(s)/=size(f6,6)) then; iflag = 705; return; end if
            if (size(x)/=size(bcoef6,1)) then; iflag = 800; return; end if
            if (size(y)/=size(bcoef6,2)) then; iflag = 801; return; end if
            if (size(z)/=size(bcoef6,3)) then; iflag = 802; return; end if
            if (size(q)/=size(bcoef6,4)) then; iflag = 803; return; end if
            if (size(r)/=size(bcoef6,5)) then; iflag = 804; return; end if
            if (size(s)/=size(bcoef6,6)) then; iflag = 805; return; end if
        end if

        status_ok = .true.
        iflag = 0

    end if

    contains

        pure subroutine check(s,n,k,x,t,ierrs,iflag,error)  !check t,x,n,k for validity

        implicit none

        character(len=1),intent(in)               :: s     !! coordinate string: 'x','y','z','q','r','s'
        integer,intent(in)              ,optional :: n     !! size of `x`
        integer,intent(in)              ,optional :: k     !! order
        real(wp),dimension(:),intent(in),optional :: x     !! abcissae vector
        real(wp),dimension(:),intent(in),optional :: t     !! knot vector size(n+k)
        integer,dimension(:),intent(in)           :: ierrs !! int error codes for n,k,x,t,size(x),size(t) checks
        integer,intent(out)                       :: iflag !! status return code
        logical,intent(out)                       :: error !! true if there was an error

        if (present(n) .and. present(k) .and. present(x) .and. present(t)) then
            call check_n('n'//s,n,x,[ierrs(1),ierrs(5)],iflag,error); if (error) return
            call check_k('k'//s,k,n,ierrs(2),iflag,error); if (error) return
            call check_x(s,n,x,ierrs(3),iflag,error); if (error) return
            if (iknot /= 0) then
                call check_t('t'//s,n,k,t,[ierrs(4),ierrs(6)],iflag,error); if (error) return
            end if
        end if

        end subroutine check

        pure subroutine check_n(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)     :: s
        integer,intent(in)              :: n
        real(wp),dimension(:),intent(in):: x     !! abcissae vector
        integer,dimension(2),intent(in) :: ierr  ![n<3 check, size(x)==n check]
        integer,intent(out)             :: iflag !! status return code
        logical,intent(out)             :: error

        if (n < 3) then
            !write(error_unit,'(A,1X,I5)') &
            !    trim(routine)//' - '//trim(s)//' is out of range: ',n
            iflag = ierr(1)
            error = .true.
        else
            if (size(x)/=n) then
                !write(error_unit,'(A,1X,I5)') &
                !    trim(routine)//' - '//trim(s)//' is not abscissa vector size'
                iflag = ierr(2)
                error = .true.
            else
                error = .false.
            end if
        end if

        end subroutine check_n

        pure subroutine check_k(s,k,n,ierr,iflag,error)

        implicit none

        character(len=*),intent(in) :: s
        integer,intent(in)          :: k
        integer,intent(in)          :: n
        integer,intent(in)          :: ierr
        integer,intent(out)         :: iflag !! status return code
        logical,intent(out)         :: error

        if ((k < 2) .or. (k >= n)) then
            !write(error_unit,'(A,1X,I5)') &
            !    trim(routine)//' - '//trim(s)//' is out of range: ',k
            iflag = ierr
            error = .true.
        else
            error = .false.
        end if

        end subroutine check_k

        pure subroutine check_x(s,n,x,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)       :: s
        integer,intent(in)                :: n
        real(wp),dimension(:),intent(in)  :: x
        integer,intent(in)                :: ierr
        integer,intent(out)               :: iflag !! status return code
        logical,intent(out)               :: error

        integer :: i

        error = .true.
        do i=2,n
            if (x(i) <= x(i-1)) then
                iflag = ierr
                !write(error_unit,'(A)') trim(routine)//' - '//trim(s)//&
                !            ' array must be strictly increasing'
                return
            end if
        end do
        error = .false.

        end subroutine check_x

        pure subroutine check_t(s,n,k,t,ierr,iflag,error)

        implicit none

        character(len=*),intent(in)       :: s
        integer,intent(in)                :: n
        integer,intent(in)                :: k
        real(wp),dimension(:),intent(in)  :: t
        integer,dimension(2),intent(in)   :: ierr  !! [non-decreasing check, size check]
        integer,intent(out)               :: iflag !! status return code
        logical,intent(out)               :: error

        integer :: i

        error = .true.

        if (size(t)/=(n+k)) then
            !write(error_unit,'(A)') trim(routine)//' - '//trim(s)//&
            !            ' array is not the correct size'
            iflag = ierr(2)
            return
        end if

        do i=2,n + k
            if (t(i) < t(i-1))  then
                iflag = ierr(1)
                !write(error_unit,'(A)') trim(routine)//' - '//trim(s)//&
                !            ' array must be non-decreasing'
                return
            end if
        end do
        error = .false.

        end subroutine check_t

    end subroutine check_inputs
!*****************************************************************************************

!*****************************************************************************************
!> dbknot chooses a knot sequence for interpolation of order k at the
!  data points x(i), i=1,..,n.  the n+k knots are placed in the array
!  t.  k knots are placed at each endpoint and not-a-knot end
!  conditions are used.  the remaining knots are placed at data points
!  if n is even and between data points if n is odd.  the rightmost
!  knot is shifted slightly to the right to insure proper interpolation
!  at x(n) (see page 350 of the reference).
!
!# History
!
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    pure subroutine dbknot(x,n,k,t)

    implicit none

    integer,intent(in)                 :: n
    integer,intent(in)                 :: k
    real(wp),dimension(n),intent(in)   :: x
    real(wp),dimension(:),intent(out)  :: t

    integer  :: i, j, ipj, npj, ip1, jstrt
    real(wp) :: rnot

    !put k knots at each endpoint
    !(shift right endpoints slightly -- see pg 350 of reference)
    rnot = x(n) + 0.1_wp*( x(n)-x(n-1) )
    do j=1,k
        t(j)   = x(1)
        npj    = n + j
        t(npj) = rnot
    end do

    !distribute remaining knots

    if (mod(k,2) == 1)  then

        !case of odd k --  knots between data points

        i = (k-1)/2 - k
        ip1 = i + 1
        jstrt = k + 1
        do j=jstrt,n
            ipj = i + j
            t(j) = 0.5_wp*( x(ipj) + x(ipj+1) )
        end do

    else

        !case of even k --  knots at data points

        i = (k/2) - k
        jstrt = k+1
        do j=jstrt,n
            ipj = i + j
            t(j) = x(ipj)
        end do

    end if

    end subroutine dbknot
!*****************************************************************************************

!*****************************************************************************************
!> dbtpcf computes b-spline interpolation coefficients for nf sets
!  of data stored in the columns of the array fcn. the b-spline
!  coefficients are stored in the rows of bcoef however.
!  each interpolation is based on the n abcissa stored in the
!  array x, and the n+k knots stored in the array t. the order
!  of each interpolation is k.
!
!# History
!
!  * Jacob Williams, 2/24/2015 : Refactored this routine.

    pure subroutine dbtpcf(x,n,fcn,ldf,nf,t,k,bcoef,work,iflag)

    integer,intent(in)                    :: n
    integer,intent(in)                    :: nf
    integer,intent(in)                    :: ldf
    integer,intent(in)                    :: k
    real(wp),dimension(n),intent(in)      :: x
    real(wp),dimension(ldf,nf),intent(in) :: fcn
    real(wp),dimension(*),intent(in)      :: t
    real(wp),dimension(nf,n),intent(out)  :: bcoef
    real(wp),dimension(*),intent(out)     :: work   !! work array of size >= `2*k*(n+1)`
    integer,intent(out)                   :: iflag  !!   0: no errors
                                                    !! 301: n should be >0

    integer :: i, j, m1, m2, iq, iw

    ! check for null input

    if (nf > 0)  then

        ! partition work array
        m1 = k - 1
        m2 = m1 + k
        iq = 1 + n
        iw = iq + m2*n+1

        ! compute b-spline coefficients

        ! first data set

        call dbintk(x,fcn,t,n,k,work,work(iq),work(iw),iflag)
        if (iflag == 0) then
            do i=1,n
                bcoef(1,i) = work(i)
            end do

            !  all remaining data sets by back-substitution

            if (nf == 1)  return
            do j=2,nf
                do i=1,n
                    work(i) = fcn(i,j)
                end do
                call dbnslv(work(iq),m2,n,m1,m1,work)
                do i=1,n
                    bcoef(j,i) = work(i)
                end do
            end do
        end if

    else
        !write(error_unit,'(A)') 'dbtpcf - n should be >0'
        iflag = 301
    end if

    end subroutine dbtpcf
!*****************************************************************************************

!*****************************************************************************************
!> dbintk produces the b-spline coefficients, bcoef, of the
!  b-spline of order k with knots t(i), i=1,...,n+k, which
!  takes on the value y(i) at x(i), i=1,...,n.  the spline or
!  any of its derivatives can be evaluated by calls to [[dbvalu]].
!
!  the i-th equation of the linear system a*bcoef = b for the
!  coefficients of the interpolant enforces interpolation at
!  x(i), i=1,...,n.  hence, b(i) = y(i), for all i, and a is
!  a band matrix with 2k-1 bands if a is invertible.  the matrix
!  a is generated row by row and stored, diagonal by diagonal,
!  in the rows of q, with the main diagonal going into row k.
!  the banded system is then solved by a call to dbnfac (which
!  constructs the triangular factorization for a and stores it
!  again in q), followed by a call to dbnslv (which then
!  obtains the solution bcoef by substitution).  dbnfac does no
!  pivoting, since the total positivity of the matrix a makes
!  this unnecessary.  the linear system to be solved is
!  (theoretically) invertible if and only if
!          t(i) < x(i) < t(i+k),        for all i.
!  equality is permitted on the left for i=1 and on the right
!  for i=n when k knots are used at x(1) or x(n).  otherwise,
!  violation of this condition is certain to lead to an error.
!
!# Error conditions
!
!  * improper input
!  * singular system of equations
!
!# History
!
!  * splint written by carl de boor [5]
!  * dbintk author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations. (jec)
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbintk(x,y,t,n,k,bcoef,q,work,iflag)

    implicit none

    integer,intent(in)                :: n      !! number of data points, n >= k
    real(wp),dimension(n),intent(in)  :: x      !! vector of length n containing data point abscissa
                                                !! in strictly increasing order.
    real(wp),dimension(n),intent(in)  :: y      !! corresponding vector of length n containing data
                                                !! point ordinates.
    real(wp),dimension(*),intent(in)  :: t      !! knot vector of length n+k
                                                !! since t(1),..,t(k) <= x(1) and t(n+1),..,t(n+k)
                                                !! >= x(n), this leaves only n-k knots (not
                                                !! necessarily x(i) values) interior to (x(1),x(n))
    integer,intent(in)                :: k      !! order of the spline, k >= 1
    real(wp),dimension(n),intent(out) :: bcoef  !! a vector of length n containing the b-spline coefficients
    real(wp),dimension(*),intent(out) :: q      !! a work vector of length (2*k-1)*n, containing
                                                !! the triangular factorization of the coefficient
                                                !! matrix of the linear system being solved.  the
                                                !! coefficients for the interpolant of an
                                                !! additional data set (x(i),yy(i)), i=1,...,n
                                                !! with the same abscissa can be obtained by loading
                                                !! yy into bcoef and then executing
                                                !! call dbnslv(q,2k-1,n,k-1,k-1,bcoef)
    real(wp),dimension(*),intent(out) :: work   !! work vector of length 2*k
    integer,intent(out)               :: iflag  !!   0: no errors.
                                                !! 100: k does not satisfy k>=1.
                                                !! 101: n does not satisfy n>=k.
                                                !! 102: x(i) does not satisfy x(i)<x(i+1) for some i.
                                                !! 103: some abscissa was not in the support of the.
                                                !! corresponding basis function and the system is singular.
                                                !! 104: the system of solver detects a singular system.
                                                !! although the theoretical conditions for a solution were satisfied.

    integer :: iwork, i, ilp1mx, j, jj, km1, kpkm2, left,lenq, np1
    real(wp) :: xi
    logical :: found

    if (k<1) then
        !write(error_unit,'(A)') 'dbintk - k does not satisfy k>=1'
        iflag = 100
        return
    end if

    if (n<k) then
        !write(error_unit,'(A)') 'dbintk - n does not satisfy n>=k'
        iflag = 101
        return
    end if

    jj = n - 1
    if (jj/=0) then
        do i=1,jj
            if (x(i)>=x(i+1)) then
                !write(error_unit,'(A)') 'dbintk - x(i) does not satisfy x(i)<x(i+1) for some i'
                iflag = 102
                return
            end if
        end do
    end if

    np1 = n + 1
    km1 = k - 1
    kpkm2 = 2*km1
    left = k
    ! zero out all entries of q
    lenq = n*(k+km1)
    do i=1,lenq
        q(i) = 0.0_wp
    end do

    ! loop over i to construct the n interpolation equations
    do i=1,n

        xi = x(i)
        ilp1mx = min(i+k,np1)
        ! find left in the closed interval (i,i+k-1) such that
        !         t(left) <= x(i) < t(left+1)
        ! matrix is singular if this is not possible
        left = max(left,i)
        if (xi<t(left)) then
            !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
            !             ' corresponding basis function and the system is singular'
            iflag = 103
            return
        end if
        found = .false.
        do
            found = (xi<t(left+1))
            if (found) exit
            left = left + 1
            if (left>=ilp1mx) exit
        end do
        if (.not. found) then
            left = left - 1
            if (xi>t(left+1)) then
                !write(error_unit,'(A)') 'dbintk - some abscissa was not in the support of the'//&
                !             ' corresponding basis function and the system is singular'
                iflag = 103
                return
            end if
        end if
        ! the i-th equation enforces interpolation at xi, hence
        ! a(i,j) = b(j,k,t)(xi), all j. only the  k  entries with  j =
        ! left-k+1,...,left actually might be nonzero. these  k  numbers
        ! are returned, in  bcoef (used for temp.storage here), by the
        ! following
        call dbspvn(t, k, k, 1, xi, left, bcoef, work, iwork, iflag)
        if (iflag/=0) return

        ! we therefore want  bcoef(j) = b(left-k+j)(xi) to go into
        ! a(i,left-k+j), i.e., into  q(i-(left+j)+2*k,(left+j)-k) since
        ! a(i+j,j)  is to go into  q(i+k,j), all i,j,  if we consider  q
        ! as a two-dim. array , with  2*k-1  rows (see comments in
        ! dbnfac). in the present program, we treat  q  as an equivalent
        ! one-dimensional array (because of fortran restrictions on
        ! dimension statements) . we therefore want  bcoef(j) to go into
        ! entry
        !     i -(left+j) + 2*k + ((left+j) - k-1)*(2*k-1)
        !            = i-left+1 + (left -k)*(2*k-1) + (2*k-2)*j
        ! of q.
        jj = i - left + 1 + (left-k)*(k+km1)
        do j=1,k
            jj = jj + kpkm2
            q(jj) = bcoef(j)
        end do

    end do

    ! obtain factorization of a, stored again in q.
    call dbnfac(q, k+km1, n, km1, km1, iflag)

    if (iflag==1) then !success
        ! solve  a*bcoef = y  by backsubstitution
        do i=1,n
            bcoef(i) = y(i)
        end do
        call dbnslv(q, k+km1, n, km1, km1, bcoef)
        iflag = 0
    else  !failure
        !write(error_unit,'(A)') 'dbintk - the system of solver detects a singular system'//&
        !             ' although the theoretical conditions for a solution were satisfied'
        iflag = 104
    end if

    end subroutine dbintk
!*****************************************************************************************

!*****************************************************************************************
!> Returns in w the LU-factorization (without pivoting) of the banded
!  matrix a of order nrow with (nbandl + 1 + nbandu) bands or diagonals
!  in the work array w .
!
!  gauss elimination without pivoting is used. the routine is
!  intended for use with matrices a which do not require row inter-
!  changes during factorization, especially for the totally
!  positive matrices which occur in spline calculations.
!  the routine should not be used for an arbitrary banded matrix.
!
!# Work array
!
! **Input**
!
!        w array of size (nroww,nrow) contains the interesting
!        part of a banded matrix  a , with the diagonals or bands of  a
!        stored in the rows of  w , while columns of  a  correspond to
!        columns of  w . this is the storage mode used in  linpack  and
!        results in efficient innermost loops.
!           explicitly,  a  has  nbandl  bands below the diagonal
!                            +     1     (main) diagonal
!                            +   nbandu  bands above the diagonal
!        and thus, with    middle = nbandu + 1,
!          a(i+j,j)  is in  w(i+middle,j)  for i=-nbandu,...,nbandl
!                                              j=1,...,nrow .
!        for example, the interesting entries of a (1,2)-banded matrix
!        of order  9  would appear in the first  1+1+2 = 4  rows of  w
!        as follows.
!                          13 24 35 46 57 68 79
!                       12 23 34 45 56 67 78 89
!                    11 22 33 44 55 66 77 88 99
!                    21 32 43 54 65 76 87 98
!
!        all other entries of  w  not identified in this way with an en-
!        try of  a  are never referenced .
!
! **Output**
!
!  * if  iflag = 1, then
!        w contains the lu-factorization of  a  into a unit lower triangu-
!        lar matrix  l  and an upper triangular matrix  u (both banded)
!        and stored in customary fashion over the corresponding entries
!        of  a . this makes it possible to solve any particular linear
!        system  a*x = b  for  x  by a
!              call dbnslv ( w, nroww, nrow, nbandl, nbandu, b )
!        with the solution x  contained in  b  on return .
!  * if  iflag = 2, then
!        one of  nrow-1, nbandl,nbandu failed to be nonnegative, or else
!        one of the potential pivots was found to be zero indicating
!        that  a  does not have an lu-factorization. this implies that
!        a  is singular in case it is totally positive .
!
!# History
!
!  * banfac written by carl de boor [5]
!  * dbnfac from CMLIB [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbnfac(w,nroww,nrow,nbandl,nbandu,iflag)

    integer,intent(in) :: nroww   !! row dimension of the work array w. must be >= nbandl + 1 + nbandu.
    integer,intent(in) :: nrow    !! matrix order
    integer,intent(in) :: nbandl  !! number of bands of a below the main diagonal
    integer,intent(in) :: nbandu  !! number of bands of a above the main diagonal
    integer,intent(out) :: iflag  !! indicating success(=1) or failure (=2)
    real(wp),dimension(nroww,nrow),intent(inout) :: w  !! work array. See header for details.

    integer :: i, ipk, j, jmax, k, kmax, middle, midmk, nrowm1
    real(wp) :: factor, pivot

    iflag = 1
    middle = nbandu + 1   ! w(middle,.) contains the main diagonal of a.
    nrowm1 = nrow - 1

    if (nrowm1 < 0) then
        iflag = 2
        return
    elseif (nrowm1 == 0) then
        if (w(middle,nrow)==0.0_wp) iflag = 2
        return
    end if

    if (nbandl<=0) then
        ! a is upper triangular. check that diagonal is nonzero .
        do i=1,nrowm1
            if (w(middle,i)==0.0_wp) then
                iflag = 2
                return
            end if
        end do
        if (w(middle,nrow)==0.0_wp) iflag = 2
        return
    end if

    if (nbandu<=0) then
        ! a is lower triangular. check that diagonal is nonzero and
        ! divide each column by its diagonal.
        do i=1,nrowm1
            pivot = w(middle,i)
            if (pivot==0.0_wp) then
                iflag = 2
                return
            end if
            jmax = min(nbandl,nrow-i)
            do j=1,jmax
                w(middle+j,i) = w(middle+j,i)/pivot
            end do
        end do
        return
    end if

    ! a is not just a triangular matrix. construct lu factorization
    do i=1,nrowm1
        ! w(middle,i)  is pivot for i-th step .
        pivot = w(middle,i)
        if (pivot==0.0_wp) then
            iflag = 2
            return
        end if
        ! jmax is the number of (nonzero) entries in column i
        ! below the diagonal.
        jmax = min(nbandl,nrow-i)
        ! divide each entry in column i below diagonal by pivot.
        do j=1,jmax
            w(middle+j,i) = w(middle+j,i)/pivot
        end do
        ! kmax is the number of (nonzero) entries in row i to
        ! the right of the diagonal.
        kmax = min(nbandu,nrow-i)
        ! subtract a(i,i+k)*(i-th column) from (i+k)-th column
        ! (below row i).
        do k=1,kmax
            ipk = i + k
            midmk = middle - k
            factor = w(midmk,ipk)
            do j=1,jmax
                w(midmk+j,ipk) = w(midmk+j,ipk) - w(middle+j,i)*factor
            end do
        end do
    end do

    ! check the last diagonal entry.
    if (w(middle,nrow)==0.0_wp) iflag = 2

    end subroutine dbnfac

!*****************************************************************************************
!> Companion routine to [[dbnfac]]. it returns the solution x of the
!  linear system a*x = b in place of b, given the lu-factorization
!  for a in the work array w from dbnfac.
!
!  (with \( a = l*u \), as stored in w), the unit lower triangular system
!  \( l(u*x) = b \) is solved for \( y = u*x \), and y stored in b. then the
!  upper triangular system \(u*x = y \) is solved for x. the calculations
!  are so arranged that the innermost loops stay within columns.
!
!# History
!
!  * banslv written by carl de boor [5]
!  * dbnslv from SLATEC library [1]
!  * Jacob Williams, 5/10/2015 : converted to free-form Fortran.

    pure subroutine dbnslv(w,nroww,nrow,nbandl,nbandu,b)

    integer,intent(in) :: nroww   !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    integer,intent(in) :: nrow    !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    integer,intent(in) :: nbandl  !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    integer,intent(in) :: nbandu  !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    real(wp),dimension(nroww,nrow),intent(in) :: w    !! describes the lu-factorization of a banded matrix a of order nrow as constructed in [[dbnfac]].
    real(wp),dimension(nrow),intent(inout) :: b  !! **in**: right side of the system to be solved
                                                 !! **out**: the solution x, of order nrow

    integer :: i, j, jmax, middle, nrowm1

    middle = nbandu + 1
    if (nrow/=1) then

        nrowm1 = nrow - 1
        if (nbandl/=0) then

            ! forward pass
            ! for i=1,2,...,nrow-1, subtract right side(i)*(i-th column of l)
            !                       from right side (below i-th row).
            do i=1,nrowm1
                jmax = min(nbandl,nrow-i)
                do j=1,jmax
                    b(i+j) = b(i+j) - b(i)*w(middle+j,i)
                end do
            end do

        end if

        ! backward pass
        ! for i=nrow,nrow-1,...,1, divide right side(i) by i-th diagonal
        !                          entry of u, then subtract right side(i)*(i-th column
        !                          of u) from right side (above i-th row).
        if (nbandu<=0) then
            ! a is lower triangular.
            do i=1,nrow
                b(i) = b(i)/w(1,i)
            end do
            return
        end if

        i = nrow
        do
            b(i) = b(i)/w(middle,i)
            jmax = min(nbandu,i-1)
            do j=1,jmax
                b(i-j) = b(i-j) - b(i)*w(middle-j,i)
            end do
            i = i - 1
            if (i<=1) exit
        end do

    end if

    b(1) = b(1)/w(middle,1)

    end subroutine dbnslv
!*****************************************************************************************

!*****************************************************************************************
!> Calculates the value of all (possibly) nonzero basis
!  functions at x of order max(jhigh,(j+1)*(index-1)), where t(k)
!  <= x <= t(n+1) and j=iwork is set inside the routine on
!  the first call when index=1.  ileft is such that t(ileft) <=
!  x < t(ileft+1).  a call to dintrv(t,n+1,x,ilo,ileft,mflag)
!  produces the proper ileft.  dbspvn calculates using the basic
!  algorithm needed in dbspvd.  if only basis functions are
!  desired, setting jhigh=k and index=1 can be faster than
!  calling dbspvd, but extra coding is required for derivatives
!  (index=2) and dbspvd is set up for this purpose.
!
!  left limiting values are set up as described in dbspvd.
!
!#Error Conditions
!
!  * improper input
!
!# History
!
!  * bsplvn written by carl de boor [5]
!  * dbspvn author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine dbspvn(t,jhigh,k,index,x,ileft,vnikx,work,iwork,iflag)

    implicit none

    real(wp),dimension(*),intent(in)  :: t        !! knot vector of length n+k, where
                                                  !! n = number of b-spline basis functions
                                                  !! n = sum of knot multiplicities-k
                                                  !! dimension t(ileft+jhigh)
    integer,intent(in)                :: jhigh    !! order of b-spline, 1 <= jhigh <= k
    integer,intent(in)                :: k        !! highest possible order
    integer,intent(in)                :: index    !! index = 1 gives basis functions of order jhigh
                                                  !!       = 2 denotes previous entry with work, iwork
                                                  !!         values saved for subsequent calls to
                                                  !!         dbspvn.
    real(wp),intent(in)               :: x        !! argument of basis functions, t(k) <= x <= t(n+1)
    integer,intent(in)                :: ileft    !! largest integer such that t(ileft) <= x < t(ileft+1)
    real(wp),dimension(k),intent(out) :: vnikx    !! vector of length k for spline values.
    real(wp),dimension(*),intent(out) :: work     !! a work vector of length 2*k
    integer,intent(out)               :: iwork    !! a work parameter.  both work and iwork contain
                                                  !! information necessary to continue for index = 2.
                                                  !! when index = 1 exclusively, these are scratch
                                                  !! variables and can be used for other purposes.
    integer,intent(out)               :: iflag    !!   0: no errors
                                                  !! 201: k does not satisfy k>=1
                                                  !! 202: jhigh does not satisfy 1<=jhigh<=k
                                                  !! 203: index is not 1 or 2
                                                  !! 204: x does not satisfy t(ileft)<=x<=t(ileft+1)

    integer :: imjp1, ipj, jp1, jp1ml, l
    real(wp) :: vm, vmprev

    ! content of j, deltam, deltap is expected unchanged between calls.
    ! work(i) = deltap(i),
    ! work(k+i) = deltam(i), i = 1,k

    if (k<1) then
        !write(error_unit,'(A)') 'dbspvn - k does not satisfy k>=1'
        iflag = 201
        return
    end if
    if (jhigh>k .or. jhigh<1) then
        !write(error_unit,'(A)') 'dbspvn - jhigh does not satisfy 1<=jhigh<=k'
        iflag = 202
        return
    end if
    if (index<1 .or. index>2) then
        !write(error_unit,'(A)') 'dbspvn - index is not 1 or 2'
        iflag = 203
        return
    end if
    if (x<t(ileft) .or. x>t(ileft+1)) then
        !write(error_unit,'(A)') 'dbspvn - x does not satisfy t(ileft)<=x<=t(ileft+1)'
        iflag = 204
        return
    end if

    iflag = 0

    if (index==1) then
        iwork = 1
        vnikx(1) = 1.0_wp
        if (iwork>=jhigh) return
    end if

    do
        ipj = ileft + iwork
        work(iwork) = t(ipj) - x
        imjp1 = ileft - iwork + 1
        work(k+iwork) = x - t(imjp1)
        vmprev = 0.0_wp
        jp1 = iwork + 1
        do l=1,iwork
            jp1ml = jp1 - l
            vm = vnikx(l)/(work(l)+work(k+jp1ml))
            vnikx(l) = vm*work(l) + vmprev
            vmprev = vm*work(k+jp1ml)
        end do
        vnikx(jp1) = vmprev
        iwork = jp1
        if (iwork>=jhigh) exit
    end do

    end subroutine dbspvn
!*****************************************************************************************

!*****************************************************************************************
!> Evaluates the b-representation (t,a,n,k) of a b-spline
!  at x for the function value on ideriv=0 or any of its
!  derivatives on ideriv=1,2,...,k-1.  right limiting values
!  (right derivatives) are returned except at the right end
!  point x=t(n+1) where left limiting values are computed.  the
!  spline is defined on t(k) <= x <= t(n+1).  dbvalu returns
!  a fatal error message when x is outside of this interval.
!
!  to compute left derivatives or left limiting values at a
!  knot t(i), replace n by i-1 and set x=t(i), i=k+1,n+1.
!
!#Error Conditions
!
!  * improper input
!
!# History
!
!  * bvalue written by carl de boor [5]
!  * dbvalu author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * 000330 modified array declarations.  (jec)
!  * Jacob Williams, 2/24/2015 : extensive refactoring of CMLIB routine.

    pure subroutine dbvalu(t,a,n,k,ideriv,x,inbv,work,iflag,val)

    implicit none

    real(wp),intent(out)             :: val     !! the interpolated value
    integer,intent(in)               :: n       !! number of b-spline coefficients.
                                                !! (sum of knot multiplicities-k)
    real(wp),dimension(:),intent(in) :: t       !! knot vector of length n+k
    real(wp),dimension(n),intent(in) :: a       !! b-spline coefficient vector of length n
    integer,intent(in)               :: k       !! order of the b-spline, k >= 1
    integer,intent(in)               :: ideriv  !! order of the derivative, 0 <= ideriv <= k-1.
                                                !! ideriv = 0 returns the b-spline value
    real(wp),intent(in)              :: x       !! argument, t(k) <= x <= t(n+1)
    integer,intent(inout)            :: inbv    !! an initialization parameter which must be set
                                                !! to 1 the first time dbvalu is called.
                                                !! inbv contains information for efficient process-
                                                !! ing after the initial call and inbv must not
                                                !! be changed by the user.  distinct splines require
                                                !! distinct inbv parameters.
    real(wp),dimension(:),intent(inout) :: work !! work vector of length at least 3*k
    integer,intent(out)              :: iflag   !!   0: no errors
                                                !! 401: k does not satisfy k>=1
                                                !! 402: n does not satisfy n>=k
                                                !! 403: ideriv does not satisfy 0<=ideriv<k
                                                !! 404: x is not greater than or equal to t(k)
                                                !! 405: x is not less than or equal to t(n+1)
                                                !! 406: a left limiting value cannot be obtained at t(k)

    integer :: i,iderp1,ihi,ihmkmj,ilo,imk,imkpj,ipj,&
               ip1,ip1mj,j,jj,j1,j2,kmider,kmj,km1,kpk,mflag
    real(wp) :: fkmj

    val = 0.0_wp

    if (k<1) then
        !write(error_unit,'(A)') 'dbvalu - k does not satisfy k>=1'
        iflag = 401
        return
    end if

    if (n<k) then
        !write(error_unit,'(A)') 'dbvalu - n does not satisfy n>=k'
        iflag = 402
        return
    end if

    if (ideriv<0 .or. ideriv>=k) then
        !write(error_unit,'(A)') 'dbvalu - ideriv does not satisfy 0<=ideriv<k'
        iflag = 403
        return
    end if

    kmider = k - ideriv

    ! find *i* in (k,n) such that t(i) <= x < t(i+1)
    ! (or, <= t(i+1) if t(i) < t(i+1) = t(n+1)).

    km1 = k - 1
    call dintrv(t, n+1, x, inbv, i, mflag)
    if (x<t(k)) then
        !write(error_unit,'(A)') 'dbvalu - x is not greater than or equal to t(k)'
        iflag = 404
        return
    end if

    if (mflag/=0) then

        if (x>t(i)) then
            !write(error_unit,'(A)') 'dbvalu - x is not less than or equal to t(n+1)'
            iflag = 405
            return
        end if

        do
            if (i==k) then
                !write(error_unit,'(A)') 'dbvalu - a left limiting value cannot be obtained at t(k)'
                iflag = 406
                return
            end if
            i = i - 1
            if (x/=t(i)) exit
        end do

    end if

    ! difference the coefficients *ideriv* times
    ! work(i) = aj(i), work(k+i) = dp(i), work(k+k+i) = dm(i), i=1.k

    imk = i - k
    do j=1,k
        imkpj = imk + j
        work(j) = a(imkpj)
    end do

    if (ideriv/=0) then
        do j=1,ideriv
            kmj = k - j
            fkmj = real(kmj,wp)
            do jj=1,kmj
                ihi = i + jj
                ihmkmj = ihi - kmj
                work(jj) = (work(jj+1)-work(jj))/(t(ihi)-t(ihmkmj))*fkmj
            end do
        end do
    end if

    ! compute value at *x* in (t(i),(t(i+1)) of ideriv-th derivative,
    ! given its relevant b-spline coeff. in aj(1),...,aj(k-ideriv).

    if (ideriv/=km1) then
        ip1 = i + 1
        kpk = k + k
        j1 = k + 1
        j2 = kpk + 1
        do j=1,kmider
            ipj = i + j
            work(j1) = t(ipj) - x
            ip1mj = ip1 - j
            work(j2) = x - t(ip1mj)
            j1 = j1 + 1
            j2 = j2 + 1
        end do
        iderp1 = ideriv + 1
        do j=iderp1,km1
            kmj = k - j
            ilo = kmj
            do jj=1,kmj
                work(jj) = (work(jj+1)*work(kpk+ilo)+work(jj)*&
                            work(k+jj))/(work(kpk+ilo)+work(k+jj))
                ilo = ilo - 1
            end do
        end do
    end if

    iflag = 0
    val = work(1)

    end subroutine dbvalu
!*****************************************************************************************

!*****************************************************************************************
!> Computes the largest integer ileft in 1 <= ileft <= lxt
!  such that xt(ileft) <= x where xt(*) is a subdivision of
!  the x interval.
!  precisely,
!
!```fortran
!         if            x < xt(1)   then ileft=1,   mflag=-1
!         if   xt(i) <= x < xt(i+1) then ileft=i,   mflag=0
!         if xt(lxt) <= x           then ileft=lxt, mflag=1
!```
!
!  that is, when multiplicities are present in the break point
!  to the left of x, the largest index is taken for ileft.
!
!# History
!
!  * interv written by carl de boor [5]
!  * dintrv author: amos, d. e., (snla) : date written 800901
!  * revision date 820801
!  * Jacob Williams, 2/24/2015 : updated to free-form Fortran.
!  * Jacob Williams, 2/17/2016 : additional refactoring (eliminated GOTOs).

    pure subroutine dintrv(xt,lxt,x,ilo,ileft,mflag)

    implicit none

    integer,intent(in)                 :: lxt    !! length of the `xt` vector
    real(wp),dimension(lxt),intent(in) :: xt     !! a knot or break point vector of length `lxt`
    real(wp),intent(in)                :: x      !! argument
    integer,intent(inout)              :: ilo    !! an initialization parameter which must be set
                                                 !! to 1 the first time the spline array `xt` is
                                                 !! processed by dintrv. `ilo` contains information for
                                                 !! efficient processing after the initial call and `ilo`
                                                 !! must not be changed by the user.  distinct splines
                                                 !! require distinct i`lo parameters.
    integer,intent(out)                :: ileft  !! largest integer satisfying `xt(ileft) <= x`
    integer,intent(out)                :: mflag  !! signals when `x` lies out of bounds

    integer :: ihi, istep, middle

    ihi = ilo + 1
    if ( ihi>=lxt ) then
        if ( x>=xt(lxt) ) then
            mflag = 1
            ileft = lxt
            return
        end if
        if ( lxt<=1 ) then
            mflag = -1
            ileft = 1
            return
        end if
        ilo = lxt - 1
        ihi = lxt
    endif

    if ( x>=xt(ihi) ) then

        ! now x >= xt(ilo). find upper bound
        istep = 1
        do
            ilo = ihi
            ihi = ilo + istep
            if ( ihi>=lxt ) then
                if ( x>=xt(lxt) ) then
                    mflag = 1
                    ileft = lxt
                    return
                end if
                ihi = lxt
            elseif ( x>=xt(ihi) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    else

        if ( x>=xt(ilo) ) then
            mflag = 0
            ileft = ilo
            return
        end if
        ! now x <= xt(ihi). find lower bound
        istep = 1
        do
            ihi = ilo
            ilo = ihi - istep
            if ( ilo<=1 ) then
                ilo = 1
                if ( x<xt(1) ) then
                    mflag = -1
                    ileft = 1
                    return
                end if
            elseif ( x<xt(ilo) ) then
                istep = istep*2
                cycle
            endif
            exit
        end do

    endif

    ! now xt(ilo) <= x < xt(ihi). narrow the interval
    do
        middle = (ilo+ihi)/2
        if ( middle==ilo ) then
            mflag = 0
            ileft = ilo
            return
        end if
        ! note. it is assumed that middle = ilo in case ihi = ilo+1
        if ( x<xt(middle) ) then
            ihi = middle
        else
            ilo = middle
        endif
    end do

    end subroutine dintrv
!*****************************************************************************************

!*****************************************************************************************
!>
!  Returns a message string associated with the status code.

    pure function get_status_message(iflag) result(msg)

    implicit none

    integer,intent(in)           :: iflag  !! return code from one of the routines
    character(len=:),allocatable :: msg    !! status message associated with the flag

    character(len=10) :: istr   !! for integer to string conversion
    integer           :: istat  !! for write statement

    select case (iflag)

    case(  0); msg='Successful execution'

    case(  1); msg='Error in evaluate_*d: class is not initialized'

    case(  2); msg='Error in db*ink: iknot out of range'
    case(  3); msg='Error in db*ink: nx out of range'
    case(  4); msg='Error in db*ink: kx out of range'
    case(  5); msg='Error in db*ink: x not strictly increasing'
    case(  6); msg='Error in db*ink: tx not non-decreasing'
    case(  7); msg='Error in db*ink: ny out of range'
    case(  8); msg='Error in db*ink: ky out of range'
    case(  9); msg='Error in db*ink: y not strictly increasing'
    case( 10); msg='Error in db*ink: ty not non-decreasing'
    case( 11); msg='Error in db*ink: nz out of range'
    case( 12); msg='Error in db*ink: kz out of range'
    case( 13); msg='Error in db*ink: z not strictly increasing'
    case( 14); msg='Error in db*ink: tz not non-decreasing'
    case( 15); msg='Error in db*ink: nq out of range'
    case( 16); msg='Error in db*ink: kq out of range'
    case( 17); msg='Error in db*ink: q not strictly increasing'
    case( 18); msg='Error in db*ink: tq not non-decreasing'
    case( 19); msg='Error in db*ink: nr out of range'
    case( 20); msg='Error in db*ink: kr out of range'
    case( 21); msg='Error in db*ink: r not strictly increasing'
    case( 22); msg='Error in db*ink: tr not non-decreasing'
    case( 23); msg='Error in db*ink: ns out of range'
    case( 24); msg='Error in db*ink: ks out of range'
    case( 25); msg='Error in db*ink: s not strictly increasing'
    case( 26); msg='Error in db*ink: ts not non-decreasing'
    case(700); msg='Error in db*ink: size(x) /= size(fcn,1)'
    case(701); msg='Error in db*ink: size(y) /= size(fcn,2)'
    case(702); msg='Error in db*ink: size(z) /= size(fcn,3)'
    case(703); msg='Error in db*ink: size(q) /= size(fcn,4)'
    case(704); msg='Error in db*ink: size(r) /= size(fcn,5)'
    case(705); msg='Error in db*ink: size(s) /= size(fcn,6)'
    case(706); msg='Error in db*ink: size(x) /= nx'
    case(707); msg='Error in db*ink: size(y) /= ny'
    case(708); msg='Error in db*ink: size(z) /= nz'
    case(709); msg='Error in db*ink: size(q) /= nq'
    case(710); msg='Error in db*ink: size(r) /= nr'
    case(711); msg='Error in db*ink: size(s) /= ns'
    case(712); msg='Error in db*ink: size(tx) /= nx+kx'
    case(713); msg='Error in db*ink: size(ty) /= ny+ky'
    case(714); msg='Error in db*ink: size(tz) /= nz+kz'
    case(715); msg='Error in db*ink: size(tq) /= nq+kq'
    case(716); msg='Error in db*ink: size(tr) /= nr+kr'
    case(717); msg='Error in db*ink: size(ts) /= ns+ks'
    case(800); msg='Error in db*ink: size(x) /= size(bcoef,1)'
    case(801); msg='Error in db*ink: size(y) /= size(bcoef,2)'
    case(802); msg='Error in db*ink: size(z) /= size(bcoef,3)'
    case(803); msg='Error in db*ink: size(q) /= size(bcoef,4)'
    case(804); msg='Error in db*ink: size(r) /= size(bcoef,5)'
    case(805); msg='Error in db*ink: size(s) /= size(bcoef,6)'

    case(100); msg='Error in dbintk: k does not satisfy k>=1'
    case(101); msg='Error in dbintk: n does not satisfy n>=k'
    case(102); msg='Error in dbintk: x(i) does not satisfy x(i)<x(i+1) for some i'
    case(103); msg='Error in dbintk: some abscissa was not in the support of the '//&
                    'corresponding basis function and the system is singular'
    case(104); msg='Error in dbintk: the system of solver detects a singular system '//&
                   'although the theoretical conditions for a solution were satisfied'

    case(201); msg='Error in dbspvn: k does not satisfy k>=1'
    case(202); msg='Error in dbspvn: jhigh does not satisfy 1<=jhigh<=k'
    case(203); msg='Error in dbspvn: index is not 1 or 2'
    case(204); msg='Error in dbspvn: x does not satisfy t(ileft)<=x<=t(ileft+1)'

    case(301); msg='Error in dbtpcf: n should be > 0'

    case(401); msg='Error in dbvalu: k does not satisfy k>=1'
    case(402); msg='Error in dbvalu: n does not satisfy n>=k'
    case(403); msg='Error in dbvalu: ideriv does not satisfy 0<=ideriv<k'
    case(404); msg='Error in dbvalu: x is not greater than or equal to t(k)'
    case(405); msg='Error in dbvalu: x is not less than or equal to t(n+1)'
    case(406); msg='Error in dbvalu: a left limiting value cannot be obtained at t(k)'

    case(501); msg='Error in initialize_*d_specify_knots: tx is not the correct size (kx+nx)'
    case(502); msg='Error in initialize_*d_specify_knots: ty is not the correct size (ky+ny)'
    case(503); msg='Error in initialize_*d_specify_knots: tz is not the correct size (kz+nz)'
    case(504); msg='Error in initialize_*d_specify_knots: tq is not the correct size (kq+nq)'
    case(505); msg='Error in initialize_*d_specify_knots: tr is not the correct size (kr+nr)'
    case(506); msg='Error in initialize_*d_specify_knots: ts is not the correct size (ks+ns)'

    case(601); msg='Error in db*val: x value out of bounds'
    case(602); msg='Error in db*val: y value out of bounds'
    case(603); msg='Error in db*val: z value out of bounds'
    case(604); msg='Error in db*val: q value out of bounds'
    case(605); msg='Error in db*val: r value out of bounds'
    case(606); msg='Error in db*val: s value out of bounds'

    case default
        write(istr,fmt='(I10)',iostat=istat) iflag
        msg = 'Unknown status flag: '//trim(adjustl(istr))
    end select

    end function get_status_message
!*****************************************************************************************

!*****************************************************************************************
    end module bspline
!*****************************************************************************************

module interface_helmod

use emission, only: dp, linspace

implicit none

integer, parameter  :: Nr=32
real(dp), parameter :: Jlim=1.d-22, Fmin = 0.5d0

type, public       :: InterData

    real(dp), allocatable   :: F(:), R(:), gamma(:), Jem(:), heat(:)
                            !list of parametera defining barrier
    real(dp)                :: W=4.5d0, kT=0.05d0, grid_spacing(3) ! General external parameters
    integer, allocatable    :: Nstart(:,:) 
                            !indices of starting surface surf_points
    real(dp)                :: rline(Nr), grid(3) 
                            !length of V lines and grid spacing
    real(dp)                :: TimeCur=0.d0, TimeInSet=0.d0, TimeFit=0.d0, TimeInt=0.d0 
                            !timing variables
end type InterData
    

 contains
 
function surf_points(phi) result(inds2)

real(dp), intent(in)        :: phi(:,:,:)
integer , allocatable       :: inds2(:,:), inds(:,:)

integer                     :: i,j,k,Nx,Ny,Nz,N,icount, sz(3)

sz=shape(phi)
Nx=sz(1)
Ny=sz(2)
Nz=sz(3)
N=Nx*Ny*Nz
allocate(inds(3,N))

icount=1
do k=2,Nz-1
    do j=2,Ny-1
        do i=2,Nx-1
            if (phi(i,j,k)<1.d-8 .and. ( &
                    phi(i+1,j,k)>1.d-8 .or. phi(i-1,j,k)>1.d-8 .or. &
                    phi(i,j+1,k)>1.d-8 .or. phi(i,j-1,k)>1.d-8 .or. &
                    phi(i,j,k+1)>1.d-8 .or. phi(i,j,k-1)>1.d-8)) &
                    then
                inds(:,icount)=[i,j,k]
                
                icount=icount+1
            endif
        enddo
    enddo
enddo

allocate(inds2(3,icount-1))
inds2=inds(:,1:icount-1)
deallocate(inds)

end function surf_points

subroutine J_from_phi(phi,this)
    use bspline, only: db3ink,db3val
    
    type(InterData), intent(inout)  :: this
    real(dp), intent(in)            :: phi(:,:,:)
    
    integer, parameter      :: kx=2, ky=2, kz=2, iknot=0
    integer                 :: inbvx=1,inbvy=1,inbvz=1,iloy=1,iloz=1
    integer                 :: idx=0, idy=0, idz=0, iflag
    real(dp), allocatable   :: fcn(:,:,:), tx(:), ty(:), tz(:), x(:), y(:), z(:)
    !the above are spline-related parameters
    
    integer                 :: nx,ny,nz, i, j, istart,jstart,kstart, Npoints
    real(dp), dimension(3)  :: direc, Efstart
    real(dp), dimension(Nr) :: xline, yline, zline, Vline
    
    real(dp)                :: t1, t2

    
    call cpu_time(t1)
    this%rline = linspace(0.d0,2.d0,Nr)
    
    nx = size(phi,1)
    ny = size(phi,2)
    nz = size(phi,3)
    
    this%Nstart = surf_points(phi)
    Npoints = size(this%Nstart,2)
    allocate(this%F(Npoints),this%R(Npoints),this%gamma(Npoints), &
            this%Jem(Npoints), this%heat(Npoints))
    

    allocate(x(nx),y(ny),z(nz),fcn(nx,ny,nz),tx(nx+kx),ty(ny+ky),tz(nz+kz))
    x = [(this%grid_spacing(1)*(i-1), i=1,nx)]
    y = [(this%grid_spacing(2)*(i-1), i=1,ny)]
    z = [(this%grid_spacing(3)*(i-1), i=1,nz)]
    call db3ink(x,nx,y,ny,z,nz,phi,kx,ky,kz,iknot,tx,ty,tz,fcn,iflag)
    call cpu_time(t2)
    this%TimeInSet = t2 - t1
    
    do j=1,Npoints
        istart = this%Nstart(1,j)
        jstart = this%Nstart(2,j)
        kstart = this%Nstart(3,j)
    
    !find direction of the line (same as Efield direction)
        Efstart=([phi(istart+1,jstart,kstart), &
                 phi(istart,jstart+1,kstart), &
                 phi(istart,jstart,kstart+1)]-[phi(istart-1,jstart,kstart), &
                 phi(istart,jstart-1,kstart), &
                 phi(istart,jstart,kstart-1)])/this%grid_spacing
                 
        if (norm2(Efstart)> Fmin) then
                 
            direc=Efstart/norm2(Efstart)
                !set the line of interpolation and interpolate
            xline=x(istart)+this%rline*direc(1)
            yline=y(jstart)+this%rline*direc(2)
            zline=z(kstart)+this%rline*direc(3)
                
            call cpu_time(t1)
            do i=1,Nr
                call db3val(xline(i),yline(i),zline(i),idx,idy,idz,tx,ty,tz,nx,ny,nz, &
                        kx,ky,kz,fcn,Vline(i),iflag,inbvx,inbvy,inbvz,iloy,iloz)
            enddo
            call cpu_time(t2)
            this%TimeInt = this%TimeInt + t2 - t1
                
            call cpu_time(t1)
            call fitpot(this%rline,Vline,this%F(j),this%R(j),this%gamma(j))
            call cpu_time(t2)
            this%TimeFit = this%TimeFit + t2 - t1
        else
            this%F(j) = norm2(Efstart)
            this%R(j) = 1.d4
            this%gamma(j) = 1.d0
        endif
            
        call cpu_time(t1)
        call current(this,j)
        call cpu_time(t2)
        this%TimeCur = this%TimeCur + t2 - t1
        

    enddo
end subroutine J_from_phi

subroutine current(this,i)
!calls emission module and calculates emission for the ith 
    use emission, only: EmissionData, cur_dens, print_data
    
    type(InterData), intent(inout)  :: this
    integer, intent(in)             :: i
    
    type(EmissionData)               :: emit

    emit%F = this%F(i)
    emit%W = this%W
    emit%R = abs(this%R(i))
    emit%gamma = this%gamma(i)
    emit%kT = this%kT
    emit%full = .false.
    
    call cur_dens(emit)
    
    if (emit%Jem > Jlim) then
        emit%full = .true.
        call cur_dens(emit)
    endif
    
    if (isnan(emit%Jem) .or. emit%Jem > 1.d201) then
        call print_data(emit)
    endif
    
    this%Jem(i) = emit%Jem
    this%heat(i) = emit%heat
    
end subroutine current

subroutine fitpot(x,V,F,R,gamma)
    !Fits V(x) data to F,R,gamma standard potential using L-V
    !minimization module
    use Levenberg_Marquardt, only: nlinfit
    
    real(dp), intent(in)    ::x(:),V(:)
    real(dp), intent(out)   ::F,R,gamma
    real(dp)                ::p(3),F2,Fend,var
    
    integer                 :: Nstart=3,i
  
    
    p(1) = (V(Nstart)-V(Nstart-1))/(x(Nstart)-x(Nstart-1))
    F2 = (V(Nstart+1)-V(Nstart))/(x(Nstart+1)-x(Nstart))
    Fend = (V(size(V))-V(size(V)-1))/(x(size(x))-x(size(x)-1))
    p(2) = abs(2.d0/((F2-p(1))/(x(Nstart)-x(Nstart-1))))
    p(3) = p(1)/Fend
!     print *, 'F,R,gamma first guess =', p
    var=nlinfit(fun,x,V,p)
    F=p(1)
    R=p(2)
    gamma=p(3)

    contains

    pure function fun(x,p) result(y)
        real(dp), intent(in) :: x
        real(dp), intent(in) :: p(:)
        real(dp)             :: y
        y=(p(1)*p(2)*x*(p(3)-1.d0)+p(1)*x**2) / (p(3)*x+p(2)*(p(3)-1.d0))
    end function fun
end subroutine fitpot

end module interface_helmod
