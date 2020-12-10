function  [displ_sp1,displddot_sp1,displdot_sp1,eps_xx_0,eps_xx_z,eps_xx_y,gam_xz_0,gam_xz_y,gam_xy_0,gam_xy_z,FReaction_sp1] = structuralDynamicsTransient(rotationalEffects,...
      gravityOn,...
      nlOn,...
      maxNumLoadSteps,...
      airDensity,...
      RayleighAlpha,...
      RayleighBeta,...
      elementOrder,...
      delta_t,...
      platformTurbineConnectionNodeNumber,...
      tolerance,...
      maxIterations,...
      numpBC,...
      numEl,...
      Omega,...
      OmegaDot,...
      mel,...
      moiel,...
      xmel,...
      Fexternal,...
      Fdof,...
      CN2H,...
      rbData,...
      jointTransform,...
      conn,...
      joint,...
      pBC,...
      x,...
      y,...
      z,...
      ac,...
      twist,...
      rhoA,...
      EIyy,...
      EIzz,...
      GJ,...
      EA,...
      rhoIyy,...
      rhoIzz,...
      rhoJ,...
      zcm,...
      ycm,...
      a,...
      EIyz,...
      alpha1,...
      alpha2,...
      alpha3,...
      alpha4,...
      alpha5,...
      alpha6,...
      rhoIyz,...
      b,...
      a0,...
      aeroCenterOffset,...
      elLen,...
      psi,...
      theta,...
      roll,...
      displ_s,...
      displdot_s,...
      displddot_s,...
      concLoadnodeNum,...
      concLoaddof,...
      concLoadval,...
      concStiffnodeNum,...
      concStiffdof,...
      concStiffval,...
      concMassnodeNum,...
      concMassdof,...
      concMassval,...
      StiffGennodeNum,...
      StiffGendof1,...
      StiffGendof2,...
      StiffGenval,...
      MassGennodeNum,...
      MassGendof1,...
      MassGendof2,...
      MassGenval,...
      DampGennodeNum,...
      DampGendof1,...
      DampGendof2,...
      DampGenval,...
      K11,...
      K12,...
      K13,...
      K14,...
      K15,...
      K16,...
      K22,...
      K23,...
      K24,...
      K25,...
      K26,...
      K33,...
      K34,...
      K35,...
      K36,...
      K44,...
      K45,...
      K46,...
      K55,...
      K56,...
      K66,...
      M11,...
      M15,...
      M16,...
      M22,...
      M24,...
      M33,...
      M34,...
      M44,...
      M55,...
      M56,...
      M66,...
      S11,...
      S12,...
      S13,...
      S15,...
      S16,...
      S22,...
      S23,...
      S25,...
      S26,...
      S33,...
      S35,...
      S36,...
      S55,...
      S56,...
      S66,...
      S14_1,...
      S14_2,...
      S24_1,...
      S24_2,...
      S34_1,...
      S34_2,...
      S45_1,...
      S45_2,...
      S46_1,...
      S46_2,...
      S44_1,...
      S44_2,...
      S44_3,...
      C12,...
      C13,...
      C23,...
      C24,...
      C25,...
      C26,...
      C34,...
      C35,...
      C36,...
      C14_1,...
      C14_2,...
      C45_1,...
      C45_2,...
      C46_1,...
      C46_2)

      % if rotationalEffects
      %       fprintf('roteffects true\n');
      % end
      % if gravityOn
      %       fprintf('gravity true\n');
      % end
      % if nlOn
      %       fprintf('nlOn true\n');
      % end

      % fprintf('maxNumLoadSteps: ');
      % fprintf('%1.1e\n',maxNumLoadSteps);
      % fprintf('airDensity: ');
      % fprintf('%1.1e\n',airDensity);
      % fprintf('RayleighAlpha: ');
      % fprintf('%1.1e\n',RayleighAlpha);
      % fprintf('RayleighBeta: ');
      % fprintf('%1.1e\n',RayleighBeta);
      % fprintf('elementOrder: ');
      % fprintf('%1.1e\n',elementOrder);
      % fprintf('delta_t: ');
      % fprintf('%1.1e\n',delta_t);
      % fprintf('platformTurbineConnectionNodeNumber: ');
      % fprintf('%1.1e\n',platformTurbineConnectionNodeNumber);
      % fprintf('tolerance: ');
      % fprintf('%1.1e\n',tolerance);
      % fprintf('maxIterations: ');
      % fprintf('%1.1e\n',maxIterations);
      % fprintf('numpBC: ');
      % fprintf('%1.1e\n',numpBC);
      % fprintf('numEl: ');
      % fprintf('%1.1e\n',numEl);
      % fprintf('Omega: ');
      % fprintf('%1.1e\n',Omega);
      % fprintf('OmegaDot: ');
      % fprintf('%1.1e\n',OmegaDot);
      % fprintf('mel: ');
      % fprintf('%1.1e\n',mel(1));
      % fprintf('moiel: ');
      % fprintf('%1.1e\n',moiel(1));
      % fprintf('xmel: ');
      % fprintf('%1.1e\n',xmel(1));
      % fprintf('Fexternal: ');
      % fprintf('%1.1e\n',Fexternal(1));
      % fprintf('Fdof: ');
      % fprintf('%1.1e\n',Fdof(1));
      % fprintf('CN2H: ');
      % fprintf('%1.1e\n',CN2H(1));
      % fprintf('rbData: ');
      % fprintf('%1.1e\n',rbData(1));
      % fprintf('jointTransform: ');
      % fprintf('%1.1e\n',jointTransform(1));
      % fprintf('conn: ');
      % fprintf('%1.1e\n',conn(1));
      % fprintf('joint: ');
      % fprintf('%1.1e\n',joint(1));
      % fprintf('pBC: ');
      % fprintf('%1.1e\n',pBC(1));
      % fprintf('x: ');
      % fprintf('%1.1e\n',x(1));
      % fprintf('y: ');
      % fprintf('%1.1e\n',y(1));
      % fprintf('z: ');
      % fprintf('%1.1e\n',z(1));
      % fprintf('ac: ');
      % fprintf('%1.1e\n',ac(1));
      % fprintf('twist: ');
      % fprintf('%1.1e\n',twist(1));
      % fprintf('rhoA: ');
      % fprintf('%1.1e\n',rhoA(1));
      % fprintf('EIyy: ');
      % fprintf('%1.1e\n',EIyy(1));
      % fprintf('EIzz: ');
      % fprintf('%1.1e\n',EIzz(1));
      % fprintf('GJ: ');
      % fprintf('%1.1e\n',GJ(1));
      % fprintf('EA: ');
      % fprintf('%1.1e\n',EA(1));
      % fprintf('rhoIyy: ');
      % fprintf('%1.1e\n',rhoIyy(1));
      % fprintf('rhoIzz: ');
      % fprintf('%1.1e\n',rhoIzz(1));
      % fprintf('rhoJ: ');
      % fprintf('%1.1e\n',rhoJ(1));
      % fprintf('zcm: ');
      % fprintf('%1.1e\n',zcm(1));
      % fprintf('ycm: ');
      % fprintf('%1.1e\n',ycm(1));
      % fprintf('a: ');
      % fprintf('%1.1e\n',a(1));
      % fprintf('EIyz: ');
      % fprintf('%1.1e\n',EIyz(1));
      % fprintf('alpha1: ');
      % fprintf('%1.1e\n',alpha1(1));
      % fprintf('alpha2: ');
      % fprintf('%1.1e\n',alpha2(1));
      % fprintf('alpha3: ');
      % fprintf('%1.1e\n',alpha3(1));
      % fprintf('alpha4: ');
      % fprintf('%1.1e\n',alpha4(1));
      % fprintf('alpha5: ');
      % fprintf('%1.1e\n',alpha5(1));
      % fprintf('alpha6: ');
      % fprintf('%1.1e\n',alpha6(1));
      % fprintf('rhoIyz: ');
      % fprintf('%1.1e\n',rhoIyz(1));
      % fprintf('b: ');
      % fprintf('%1.1e\n',b(1));
      % fprintf('a0: ');
      % fprintf('%1.1e\n',a0(1));
      % fprintf('aeroCenterOffset: ');
      % fprintf('%1.1e\n',aeroCenterOffset(1));
      % fprintf('elLen: ');
      % fprintf('%1.1e\n',elLen(1));
      % fprintf('psi: ');
      % fprintf('%1.1e\n',psi(1));
      % fprintf('theta: ');
      % fprintf('%1.1e\n',theta(1));
      % fprintf('roll: ');
      % fprintf('%1.1e\n',roll(1));
      % fprintf('displ_s: ');
      % fprintf('%1.1e\n',displ_s(1));
      % fprintf('displdot_s: ');
      % fprintf('%1.1e\n',displdot_s(1));
      % fprintf('displddot_s: ');
      % fprintf('%1.1e\n',displddot_s(1));
      % fprintf('K11: ');
      % fprintf('%1.1e\n',K11(1));
      % fprintf('K12: ');
      % fprintf('%1.1e\n',K12(1));
      % fprintf('K13: ');
      % fprintf('%1.1e\n',K13(1));
      % fprintf('K14: ');
      % fprintf('%1.1e\n',K14(1));
      % fprintf('K15: ');
      % fprintf('%1.1e\n',K15(1));
      % fprintf('K16: ');
      % fprintf('%1.1e\n',K16(1));
      % fprintf('K22: ');
      % fprintf('%1.1e\n',K22(1));
      % fprintf('K23: ');
      % fprintf('%1.1e\n',K23(1));
      % fprintf('K24: ');
      % fprintf('%1.1e\n',K24(1));
      % fprintf('K25: ');
      % fprintf('%1.1e\n',K25(1));
      % fprintf('K26: ');
      % fprintf('%1.1e\n',K26(1));
      % fprintf('K33: ');
      % fprintf('%1.1e\n',K33(1));
      % fprintf('K34: ');
      % fprintf('%1.1e\n',K34(1));
      % fprintf('K35: ');
      % fprintf('%1.1e\n',K35(1));
      % fprintf('K36: ');
      % fprintf('%1.1e\n',K36(1));
      % fprintf('K44: ');
      % fprintf('%1.1e\n',K44(1));
      % fprintf('K45: ');
      % fprintf('%1.1e\n',K45(1));
      % fprintf('K46: ');
      % fprintf('%1.1e\n',K46(1));
      % fprintf('K55: ');
      % fprintf('%1.1e\n',K55(1));
      % fprintf('K56: ');
      % fprintf('%1.1e\n',K56(1));
      % fprintf('K66: ');
      % fprintf('%1.1e\n',K66(1));
      % fprintf('M11: ');
      % fprintf('%1.1e\n',M11(1));
      % fprintf('M15: ');
      % fprintf('%1.1e\n',M15(1));
      % fprintf('M16: ');
      % fprintf('%1.1e\n',M16(1));
      % fprintf('M22: ');
      % fprintf('%1.1e\n',M22(1));
      % fprintf('M24: ');
      % fprintf('%1.1e\n',M24(1));
      % fprintf('M33: ');
      % fprintf('%1.1e\n',M33(1));
      % fprintf('M34: ');
      % fprintf('%1.1e\n',M34(1));
      % fprintf('M44: ');
      % fprintf('%1.1e\n',M44(1));
      % fprintf('M55: ');
      % fprintf('%1.1e\n',M55(1));
      % fprintf('M56: ');
      % fprintf('%1.1e\n',M56(1));
      % fprintf('M66: ');
      % fprintf('%1.1e\n',M66(1));
      % fprintf('S11: ');
      % fprintf('%1.1e\n',S11(1));
      % fprintf('S12: ');
      % fprintf('%1.1e\n',S12(1));
      % fprintf('S13: ');
      % fprintf('%1.1e\n',S13(1));
      % fprintf('S15: ');
      % fprintf('%1.1e\n',S15(1));
      % fprintf('S16: ');
      % fprintf('%1.1e\n',S16(1));
      % fprintf('S22: ');
      % fprintf('%1.1e\n',S22(1));
      % fprintf('S23: ');
      % fprintf('%1.1e\n',S23(1));
      % fprintf('S25: ');
      % fprintf('%1.1e\n',S25(1));
      % fprintf('S26: ');
      % fprintf('%1.1e\n',S26(1));
      % fprintf('S33: ');
      % fprintf('%1.1e\n',S33(1));
      % fprintf('S35: ');
      % fprintf('%1.1e\n',S35(1));
      % fprintf('S36: ');
      % fprintf('%1.1e\n',S36(1));
      % fprintf('S55: ');
      % fprintf('%1.1e\n',S55(1));
      % fprintf('S56: ');
      % fprintf('%1.1e\n',S56(1));
      % fprintf('S66: ');
      % fprintf('%1.1e\n',S66(1));
      % fprintf('S14_1: ');
      % fprintf('%1.1e\n',S14_1(1));
      % fprintf('S14_2: ');
      % fprintf('%1.1e\n',S14_2(1));
      % fprintf('S24_1: ');
      % fprintf('%1.1e\n',S24_1(1));
      % fprintf('S24_2: ');
      % fprintf('%1.1e\n',S24_2(1));
      % fprintf('S34_1: ');
      % fprintf('%1.1e\n',S34_1(1));
      % fprintf('S34_2: ');
      % fprintf('%1.1e\n',S34_2(1));
      % fprintf('S45_1: ');
      % fprintf('%1.1e\n',S45_1(1));
      % fprintf('S45_2: ');
      % fprintf('%1.1e\n',S45_2(1));
      % fprintf('S46_1: ');
      % fprintf('%1.1e\n',S46_1(1));
      % fprintf('S46_2: ');
      % fprintf('%1.1e\n',S46_2(1));
      % fprintf('S44_1: ');
      % fprintf('%1.1e\n',S44_1(1));
      % fprintf('S44_2: ');
      % fprintf('%1.1e\n',S44_2(1));
      % fprintf('S44_3: ');
      % fprintf('%1.1e\n',S44_3(1));
      % fprintf('C12: ');
      % fprintf('%1.1e\n',C12(1));
      % fprintf('C13: ');
      % fprintf('%1.1e\n',C13(1));
      % fprintf('C23: ');
      % fprintf('%1.1e\n',C23(1));
      % fprintf('C24: ');
      % fprintf('%1.1e\n',C24(1));
      % fprintf('C25: ');
      % fprintf('%1.1e\n',C25(1));
      % fprintf('C26: ');
      % fprintf('%1.1e\n',C26(1));
      % fprintf('C34: ');
      % fprintf('%1.1e\n',C34(1));
      % fprintf('C35: ');
      % fprintf('%1.1e\n',C35(1));
      % fprintf('C36: ');
      % fprintf('%1.1e\n',C36(1));
      % fprintf('C14_1: ');
      % fprintf('%1.1e\n',C14_1(1));
      % fprintf('C14_2: ');
      % fprintf('%1.1e\n',C14_2(1));
      % fprintf('C45_1: ');
      % fprintf('%1.1e\n',C45_1(1));
      % fprintf('C45_2: ');
      % fprintf('%1.1e\n',C45_2(1));
      % fprintf('C46_1: ');
      % fprintf('%1.1e\n',C46_1(1));
      % fprintf('C46_2: ');
      % fprintf('%1.1e\n',C46_2(1));

      model.analysisType = 'TNB';


      model.delta_t = delta_t;


      model.airDensity = airDensity;

      model.gravityOn = gravityOn;
      model.nlOn = nlOn;
      model.RayleighAlpha = RayleighAlpha;
      model.RayleighBeta = RayleighBeta;
      model.elementOrder = elementOrder;
      model.joint = joint;
      model.jointTransform = jointTransform;



      bladeData.nodeNum = 0.0;
      model.bladeData = bladeData; %Struct



      nlParams.adaptiveLoadSteppingFlag = false;
      nlParams.maxNumLoadSteps = maxNumLoadSteps;
      nlParams.tolerance = tolerance;
      nlParams.maxIterations = maxIterations;
      model.nlParams = nlParams; %Struct



      BC.numpBC = numpBC;
      BC.pBC = pBC;
      model.BC = BC; %Struct

      single_concstruct.nodeNum = 1.0;
      single_concstruct.dof = 1.0;
      single_concstruct.val = 1.0;

      concLoad = repmat(single_concstruct,1,length(concLoadnodeNum));
      for ii = 1:length(concLoadnodeNum)
            concLoad(ii).nodeNum = concLoadnodeNum(ii);
            concLoad(ii).dof = concLoaddof(ii);
            concLoad(ii).val = concLoadval(ii);
      end

      concStiff = repmat(single_concstruct,1,length(concStiffnodeNum));
      for ii = 1:length(concStiffnodeNum)
            concStiff(ii).nodeNum = concStiffnodeNum(ii);
            concStiff(ii).dof = concStiffdof(ii);
            concStiff(ii).val = concStiffval(ii);
      end

      concMass = repmat(single_concstruct,1,length(concMassnodeNum));
      for ii = 1:length(concMassnodeNum)
            concMass(ii).nodeNum = concMassnodeNum(ii);
            concMass(ii).dof = concMassdof(ii);
            concMass(ii).val = concMassval(ii);
      end

      single_genstruct.nodeNum = 1.0;
      single_genstruct.dof1 = 1.0;
      single_genstruct.dof2 = 1.0;
      single_genstruct.val = 1.0;

      concStiffGen = repmat(single_genstruct,1,length(StiffGennodeNum));
      for ii = 1:length(StiffGennodeNum)
            concStiffGen(ii).nodeNum = StiffGennodeNum(ii);
            concStiffGen(ii).dof1 = StiffGendof1(ii);
            concStiffGen(ii).dof2 = StiffGendof2(ii);
            concStiffGen(ii).val = StiffGenval(ii);
      end

      concMassGen = repmat(single_genstruct,1,length(MassGennodeNum));
      for ii = 1:length(MassGennodeNum)
            concMassGen(ii).nodeNum = MassGennodeNum(ii);
            concMassGen(ii).dof1 = MassGendof1(ii);
            concMassGen(ii).dof2 = MassGendof2(ii);
            concMassGen(ii).val = MassGenval(ii);
      end

      concDampGen = repmat(single_genstruct,1,length(DampGennodeNum));
      for ii = 1:length(DampGennodeNum)
            concDampGen(ii).nodeNum = DampGennodeNum(ii);
            concDampGen(ii).dof1 = DampGendof1(ii);
            concDampGen(ii).dof2 = DampGendof2(ii);
            concDampGen(ii).val = DampGenval(ii);
      end

      nodalTerms.concLoad = concLoad;
      nodalTerms.concStiff = concStiff;
      nodalTerms.concMass = concMass;
      nodalTerms.concStiffGen = concStiffGen;
      nodalTerms.concMassGen = concMassGen;
      nodalTerms.concDampGen = concDampGen;
      model.nodalTerms = nodalTerms; %Struct
      % Reconstruct mesh structure

      mesh.numEl = numEl;
      mesh.x = x;
      mesh.y = y;
      mesh.z = z;
      mesh.conn = conn;



      %Reconstruct el
      single_sectionPropsArray = struct('ac',zeros(1,2),...
      'twist',zeros(1,2),...
      'rhoA',zeros(1,2),...
      'EIyy',zeros(1,2),...
      'EIzz',zeros(1,2),...
      'GJ',zeros(1,2),...
      'EA',zeros(1,2),...
      'rhoIyy',zeros(1,2),...
      'rhoIzz',zeros(1,2),...
      'rhoJ',zeros(1,2),...
      'zcm',zeros(1,2),...
      'ycm',zeros(1,2),...
      'a',zeros(1,2),...
      'EIyz',zeros(1,2),...
      'alpha1',zeros(1,2),...
      'alpha2',zeros(1,2),...
      'alpha3',zeros(1,2),...
      'alpha4',zeros(1,2),...
      'alpha5',zeros(1,2),...
      'alpha6',zeros(1,2),...
      'rhoIyz',zeros(1,2),...
      'b',zeros(1,2),...
      'a0',zeros(1,2),...
      'aeroCenterOffset',zeros(1,2));

      props = repmat(single_sectionPropsArray,1,mesh.numEl);

      for jj = 1:mesh.numEl
            props(jj).ac = ac(jj,:);
            props(jj).twist = twist(jj,:);
            props(jj).rhoA = rhoA(jj,:);
            props(jj).EIyy = EIyy(jj,:);
            props(jj).EIzz = EIzz(jj,:);
            props(jj).GJ = GJ(jj,:);
            props(jj).EA = EA(jj,:);
            props(jj).rhoIyy = rhoIyy(jj,:);
            props(jj).rhoIzz = rhoIzz(jj,:);
            props(jj).rhoJ = rhoJ(jj,:);
            props(jj).zcm = zcm(jj,:);
            props(jj).ycm = ycm(jj,:);
            props(jj).a = a(jj,:);
            props(jj).EIyz = EIyz(jj,:);
            props(jj).alpha1 = alpha1(jj,:);
            props(jj).alpha2 = alpha2(jj,:);
            props(jj).alpha3 = alpha3(jj,:);
            props(jj).alpha4 = alpha4(jj,:);
            props(jj).alpha5 = alpha5(jj,:);
            props(jj).alpha6 = alpha6(jj,:);
            props(jj).rhoIyz = rhoIyz(jj,:);
            props(jj).b = b(jj,:);
            props(jj).a0 = a0(jj,:);
            props(jj).aeroCenterOffset = aeroCenterOffset(jj,:);
      end
      el.props = props; %struct



      el.elLen = elLen;
      el.psi = psi;
      el.theta = theta;
      el.roll = roll;
      el.rotationalEffects = rotationalEffects;



      % Reconstruct dispData
      dispData.displ_s = displ_s;
      dispData.displdot_s = displdot_s;
      dispData.displddot_s = displddot_s;



      % Reconstruct elStorage
      single_elStorage = struct('K11', zeros(2,2),...
      'K12', zeros(2,2),...
      'K13', zeros(2,2),...
      'K14', zeros(2,2),...
      'K15', zeros(2,2),...
      'K16', zeros(2,2),...
      'K22', zeros(2,2),...
      'K23', zeros(2,2),...
      'K24', zeros(2,2),...
      'K25', zeros(2,2),...
      'K26', zeros(2,2),...
      'K33', zeros(2,2),...
      'K34', zeros(2,2),...
      'K35', zeros(2,2),...
      'K36', zeros(2,2),...
      'K44', zeros(2,2),...
      'K45', zeros(2,2),...
      'K46', zeros(2,2),...
      'K55', zeros(2,2),...
      'K56', zeros(2,2),...
      'K66', zeros(2,2),...
      'M11', zeros(2,2),...
      'M15', zeros(2,2),...
      'M16', zeros(2,2),...
      'M22', zeros(2,2),...
      'M24', zeros(2,2),...
      'M33', zeros(2,2),...
      'M34', zeros(2,2),...
      'M44', zeros(2,2),...
      'M55', zeros(2,2),...
      'M56', zeros(2,2),...
      'M66', zeros(2,2),...
      'S11', zeros(2,2),...
      'S12', zeros(2,2),...
      'S13', zeros(2,2),...
      'S15', zeros(2,2),...
      'S16', zeros(2,2),...
      'S22', zeros(2,2),...
      'S23', zeros(2,2),...
      'S25', zeros(2,2),...
      'S26', zeros(2,2),...
      'S33', zeros(2,2),...
      'S35', zeros(2,2),...
      'S36', zeros(2,2),...
      'S55', zeros(2,2),...
      'S56', zeros(2,2),...
      'S66', zeros(2,2),...
      'S14_1', zeros(2,2),...
      'S14_2', zeros(2,2),...
      'S24_1', zeros(2,2),...
      'S24_2', zeros(2,2),...
      'S34_1', zeros(2,2),...
      'S34_2', zeros(2,2),...
      'S45_1', zeros(2,2),...
      'S45_2', zeros(2,2),...
      'S46_1', zeros(2,2),...
      'S46_2', zeros(2,2),...
      'S44_1', zeros(2,2),...
      'S44_2', zeros(2,2),...
      'S44_3', zeros(2,2),...
      'C12', zeros(2,2),...
      'C13', zeros(2,2),...
      'C23', zeros(2,2),...
      'C24', zeros(2,2),...
      'C25', zeros(2,2),...
      'C26', zeros(2,2),...
      'C34', zeros(2,2),...
      'C35', zeros(2,2),...
      'C36', zeros(2,2),...
      'C14_1', zeros(2,2),...
      'C14_2', zeros(2,2),...
      'C45_1', zeros(2,2),...
      'C45_2', zeros(2,2),...
      'C46_1', zeros(2,2),...
      'C46_2', zeros(2,2),...
      'mel', zeros(1,1),...
      'moiel', zeros(3,3),...
      'xmel', zeros(3,1));

      elStorage = repmat(single_elStorage,1,mesh.numEl);

      for jj = 1:mesh.numEl
            elStorage(jj).K11(:,:) = [K11(jj,1:2);K11(jj,3:4)];
            elStorage(jj).K12(:,:) = [K12(jj,1:2);K12(jj,3:4)];
            elStorage(jj).K13(:,:) = [K13(jj,1:2);K13(jj,3:4)];
            elStorage(jj).K14(:,:) = [K14(jj,1:2);K14(jj,3:4)];
            elStorage(jj).K15(:,:) = [K15(jj,1:2);K15(jj,3:4)];
            elStorage(jj).K16(:,:) = [K16(jj,1:2);K16(jj,3:4)];
            elStorage(jj).K22(:,:) = [K22(jj,1:2);K22(jj,3:4)];
            elStorage(jj).K23(:,:) = [K23(jj,1:2);K23(jj,3:4)];
            elStorage(jj).K24(:,:) = [K24(jj,1:2);K24(jj,3:4)];
            elStorage(jj).K25(:,:) = [K25(jj,1:2);K25(jj,3:4)];
            elStorage(jj).K26(:,:) = [K26(jj,1:2);K26(jj,3:4)];
            elStorage(jj).K33(:,:) = [K33(jj,1:2);K33(jj,3:4)];
            elStorage(jj).K34(:,:) = [K34(jj,1:2);K34(jj,3:4)];
            elStorage(jj).K35(:,:) = [K35(jj,1:2);K35(jj,3:4)];
            elStorage(jj).K36(:,:) = [K36(jj,1:2);K36(jj,3:4)];
            elStorage(jj).K44(:,:) = [K44(jj,1:2);K44(jj,3:4)];
            elStorage(jj).K45(:,:) = [K45(jj,1:2);K45(jj,3:4)];
            elStorage(jj).K46(:,:) = [K46(jj,1:2);K46(jj,3:4)];
            elStorage(jj).K55(:,:) = [K55(jj,1:2);K55(jj,3:4)];
            elStorage(jj).K56(:,:) = [K56(jj,1:2);K56(jj,3:4)];
            elStorage(jj).K66(:,:) = [K66(jj,1:2);K66(jj,3:4)];
            elStorage(jj).M11(:,:) = [M11(jj,1:2);M11(jj,3:4)];
            elStorage(jj).M15(:,:) = [M15(jj,1:2);M15(jj,3:4)];
            elStorage(jj).M16(:,:) = [M16(jj,1:2);M16(jj,3:4)];
            elStorage(jj).M22(:,:) = [M22(jj,1:2);M22(jj,3:4)];
            elStorage(jj).M24(:,:) = [M24(jj,1:2);M24(jj,3:4)];
            elStorage(jj).M33(:,:) = [M33(jj,1:2);M33(jj,3:4)];
            elStorage(jj).M34(:,:) = [M34(jj,1:2);M34(jj,3:4)];
            elStorage(jj).M44(:,:) = [M44(jj,1:2);M44(jj,3:4)];
            elStorage(jj).M55(:,:) = [M55(jj,1:2);M55(jj,3:4)];
            elStorage(jj).M56(:,:) = [M56(jj,1:2);M56(jj,3:4)];
            elStorage(jj).M66(:,:) = [M66(jj,1:2);M66(jj,3:4)];
            elStorage(jj).S11(:,:) = [S11(jj,1:2);S11(jj,3:4)];
            elStorage(jj).S12(:,:) = [S12(jj,1:2);S12(jj,3:4)];
            elStorage(jj).S13(:,:) = [S13(jj,1:2);S13(jj,3:4)];
            elStorage(jj).S15(:,:) = [S15(jj,1:2);S15(jj,3:4)];
            elStorage(jj).S16(:,:) = [S16(jj,1:2);S16(jj,3:4)];
            elStorage(jj).S22(:,:) = [S22(jj,1:2);S22(jj,3:4)];
            elStorage(jj).S23(:,:) = [S23(jj,1:2);S23(jj,3:4)];
            elStorage(jj).S25(:,:) = [S25(jj,1:2);S25(jj,3:4)];
            elStorage(jj).S26(:,:) = [S26(jj,1:2);S26(jj,3:4)];
            elStorage(jj).S33(:,:) = [S33(jj,1:2);S33(jj,3:4)];
            elStorage(jj).S35(:,:) = [S35(jj,1:2);S35(jj,3:4)];
            elStorage(jj).S36(:,:) = [S36(jj,1:2);S36(jj,3:4)];
            elStorage(jj).S55(:,:) = [S55(jj,1:2);S55(jj,3:4)];
            elStorage(jj).S56(:,:) = [S56(jj,1:2);S56(jj,3:4)];
            elStorage(jj).S66(:,:) = [S66(jj,1:2);S66(jj,3:4)];
            elStorage(jj).S14_1(:,:) = [S14_1(jj,1:2);S14_1(jj,3:4)];
            elStorage(jj).S14_2(:,:) = [S14_2(jj,1:2);S14_2(jj,3:4)];
            elStorage(jj).S24_1(:,:) = [S24_1(jj,1:2);S24_1(jj,3:4)];
            elStorage(jj).S24_2(:,:) = [S24_2(jj,1:2);S24_2(jj,3:4)];
            elStorage(jj).S34_1(:,:) = [S34_1(jj,1:2);S34_1(jj,3:4)];
            elStorage(jj).S34_2(:,:) = [S34_2(jj,1:2);S34_2(jj,3:4)];
            elStorage(jj).S45_1(:,:) = [S45_1(jj,1:2);S45_1(jj,3:4)];
            elStorage(jj).S45_2(:,:) = [S45_2(jj,1:2);S45_2(jj,3:4)];
            elStorage(jj).S46_1(:,:) = [S46_1(jj,1:2);S46_1(jj,3:4)];
            elStorage(jj).S46_2(:,:) = [S46_2(jj,1:2);S46_2(jj,3:4)];
            elStorage(jj).S44_1(:,:) = [S44_1(jj,1:2);S44_1(jj,3:4)];
            elStorage(jj).S44_2(:,:) = [S44_2(jj,1:2);S44_2(jj,3:4)];
            elStorage(jj).S44_3(:,:) = [S44_3(jj,1:2);S44_3(jj,3:4)];
            elStorage(jj).C12(:,:) = [C12(jj,1:2);C12(jj,3:4)];
            elStorage(jj).C13(:,:) = [C13(jj,1:2);C13(jj,3:4)];
            elStorage(jj).C23(:,:) = [C23(jj,1:2);C23(jj,3:4)];
            elStorage(jj).C24(:,:) = [C24(jj,1:2);C24(jj,3:4)];
            elStorage(jj).C25(:,:) = [C25(jj,1:2);C25(jj,3:4)];
            elStorage(jj).C26(:,:) = [C26(jj,1:2);C26(jj,3:4)];
            elStorage(jj).C34(:,:) = [C34(jj,1:2);C34(jj,3:4)];
            elStorage(jj).C35(:,:) = [C35(jj,1:2);C35(jj,3:4)];
            elStorage(jj).C36(:,:) = [C36(jj,1:2);C36(jj,3:4)];
            elStorage(jj).C14_1(:,:) = [C14_1(jj,1:2);C14_1(jj,3:4)];
            elStorage(jj).C14_2(:,:) = [C14_2(jj,1:2);C14_2(jj,3:4)];
            elStorage(jj).C45_1(:,:) = [C45_1(jj,1:2);C45_1(jj,3:4)];
            elStorage(jj).C45_2(:,:) = [C45_2(jj,1:2);C45_2(jj,3:4)];
            elStorage(jj).C46_1(:,:) = [C46_1(jj,1:2);C46_1(jj,3:4)];
            elStorage(jj).C46_2(:,:) = [C46_2(jj,1:2);C46_2(jj,3:4)];
            elStorage(jj).mel = mel(jj);
            elStorage(jj).moiel(:,:) = [moiel(jj,1:3);moiel(jj,4:6);moiel(jj,7:9)];
            elStorage(jj).xmel(:) = xmel(jj,:);
      end


      % **********************************************************************
      % *                   Part of the SNL OWENS Toolkit                    *
      % * Developed by Sandia National Laboratories Wind Energy Technologies *
      % *             See license.txt for disclaimer information             *
      % **********************************************************************%   [dispOut,FReaction_sp1] = structuralDynamicsTransient(model,mesh,el,...
      %                             dispData,Omega,OmegaDot,time,delta_t,...
      %                             elStorage,Fexternal,Fdof,CN2H,rbData)
      %
      %   This function performs transient structural dynamics analysis.
      %
      %   input:
      %   model      = object containing model data
      %   mesh       = object containing mesh data
      %   el         = object containing element data
      %   dispData   = object containing displacement data
      %   Omega      = rotor speed (Hz)
      %   OmegaDot   = rotor acceleratin (Hz)
      %   time       = current simulation time
      %   delta_t    = time step size
      %   elStorage  = object containing stored element data
      %   Fexternal  = vector containing external force values
      %   Fdof       = vector containing global DOF numbering associated with
      %                external force values
      %   CN2H       = transformation matrix from inertial frame to hub frame
      %   rbData     = vector containing rigid body displacement, velocity, and
      %                acceleration
      %
      %   output:
      %   dispOut       = object containing displacement data at end of time step
      %   FReaction_sp1 = vector containing reaction force at turbine base at
      %                   end of time step

      %-------- get model information -----------
      conin = conn(1,:);


      numEl = mesh.numEl;
      x = mesh.x;
      y = mesh.y;
      z = mesh.z;
      conn = mesh.conn;
      numNodes = length(x);
      elementOrder = model.elementOrder;
      BC = model.BC;

      numNodesPerEl = elementOrder + 1;
      numDOFPerNode = 6;
      totalNumDOF = numNodes * numDOFPerNode;
      % [~,numReducedDOF]=size(model.jointTransform);
      nodalTerms = model.nodalTerms;
      nodalTermsCopy = nodalTerms;
      %-----------------------------------------

      %initialize displacements, tolerance, uNorm, iteration count for nonlinear
      %iteration
      unorm = 1e6;
      tol = model.nlParams.tolerance;
      maxIterations = model.nlParams.maxIterations;
      iterationCount = 0;

      elx=zeros(numNodesPerEl,1);
      ely=zeros(numNodesPerEl,1);
      elz=zeros(numNodesPerEl,2);
      eldisp = zeros(1,numNodesPerEl*numDOFPerNode);
      eldisp_sm1 = zeros(1,numNodesPerEl*numDOFPerNode);
      eldispdot = eldisp;
      eldispddot = eldisp;
      eldispiter = eldisp;
      % if(model.nlOn)
      %      iterationType = model.nlParams.iterationType;
      % else
      %      iterationType = 'LINEAR';
      % end

      iterationType = 'DI';
      analysisType = model.analysisType;
      timeInt = struct('delta_t',0.0,...
      'a1',0.0,...
      'a2',0.0,...
      'a3',0.0,...
      'a4',0.0,...
      'a5',0.0,...
      'a6',0.0,...
      'a7',0.0,...
      'a8',0.0);
      if(strcmp(analysisType,'TNB'))
            %------ newmark integration parameters ---------
            alpha = 0.5;
            gamma = 0.5;
            beta = 0.5*gamma;

            timeInt.delta_t = delta_t;
            timeInt.a1 = alpha*delta_t;
            timeInt.a2 = (1.0-alpha)*delta_t;
            a3 = 1.0/(beta*delta_t*delta_t);
            timeInt.a3 = a3;
            timeInt.a4 = a3*delta_t;
            timeInt.a5 = 1.0/gamma-1.0;
            timeInt.a6 = alpha/(beta*delta_t);
            timeInt.a7 = alpha/beta - 1.0;
            timeInt.a8 = delta_t*(alpha/gamma-1.0);

            disp_s = dispData.displ_s;
            dispdot_s = dispData.displdot_s;
            dispddot_s = dispData.displddot_s;

            displddot_im1 = dispddot_s;
            displdot_im1 = dispdot_s;
            displ_im1 = disp_s;
      elseif(strcmp(analysisType,'TD'))
            %------ dean integration parameters -------------
            alpha = 0.25;

            timeInt.delta_t = delta_t;
            timeInt.a1 = alpha*delta_t^2;
            timeInt.a2 = (1-2*alpha)*delta_t^2;
            timeInt.a3 = delta_t/2.0;
            timeInt.a4 = delta_t*delta_t;
            disp_s = dispData.displ_s;
            %     disp_sm1 = dispData.displ_sm1;
            %-------------------------------------------
      else
            error('analysis type not supported, choose another')
      end

      %-----------------------------------------------
      % Initialize elInput, and DO NOT redundantly re-assign the memory in the
      % while and for loops below.
      elInput = struct('analysisType', analysisType,...
      'elementOrder', elementOrder,...
      'modalFlag', true,...
      'timeInt', timeInt,...
      'xloc', zeros(1,2),...
      'sectionProps', el.props(1),...
      'sweepAngle', 0.0,...
      'coneAngle', 0.0,...
      'rollAngle', 0.0,...
      'aeroSweepAngle', 0.0,...
      'firstIteration', true,...
      'concMass', zeros(4,2),...
      'concStiff', zeros(6,2),...
      'concLoad', zeros(6,2),...
      'disp', zeros(1,12),...
      'dispm1', eldisp_sm1,...
      'dispdot', eldispdot,...
      'dispddot', eldispddot,...
      'x', elx,...
      'y', ely,...
      'z', elz,...
      'accelVec', rbData(1:3),...
      'Omega', Omega,...
      'OmegaDot', OmegaDot,...
      'omegaVec', rbData(4:6),...
      'omegaDotVec', rbData(7:9),...
      'displ_iter', eldispiter,...
      'useDisp', model.nlOn,...
      'preStress', false,...
      'iterationType', iterationType,...
      'freq', 0.0,... %Is not used for this model type, but must be declared
      'aeroElasticOn', false,...
      'aeroForceOn', false,...
      'airDensity', model.airDensity,...
      'gravityOn', model.gravityOn,...
      'RayleighAlpha', model.RayleighAlpha,...
      'RayleighBeta', model.RayleighBeta,...
      'CN2H', CN2H);

      elInput.freq = 0.0; %Try to force matlab to write the C++ code so as to allocate the memory here and not inside the convergence loop.

      while(unorm>tol && iterationCount < maxIterations) %iteration loop
            %------- intitialization -----------------
            Kg = zeros(totalNumDOF,totalNumDOF); %initialize global stiffness and force vector
            Fg = zeros(totalNumDOF,1);

            nodalTerms = nodalTermsCopy;
            %-------------------------------------------

            %---- element  calculation and assembly ----------------------------------
            for i=1:numEl
                  %Calculate Ke and Fe for element i
                  index = 1;                           %initialize element data
                  elInput.xloc = [0.0 el.elLen(i)];
                  elInput.sectionProps = el.props(i);
                  elInput.sweepAngle = el.psi(i);
                  elInput.coneAngle = el.theta(i);
                  elInput.rollAngle = el.roll(i);
                  if(iterationCount == 0)
                        elInput.firstIteration = true;
                  else
                        elInput.firstIteration = false;
                  end

                  for j=1:numNodesPerEl

                        %get element cooridnates
                        elx(j) = x(conn(i,j));
                        ely(j) = y(conn(i,j));
                        elz(j) = z(conn(i,j));

                        %get element nodal displacements at s and s-1 time step
                        for k=1:numDOFPerNode
                              %                 if(strcmp(analysisType,'TD'))
                              %                     eldisp(index) = disp_s((conn(i,j)-1)*numDOFPerNode + k);
                              %                     eldisp_sm1(index) = disp_sm1((conn(i,j)-1)*numDOFPerNode + k);
                              %                     eldispiter(index) = displ_iter((conn(i,j)-1)*numDOFPerNode + k);
                              %                 end
                              if(strcmp(analysisType,'TNB'))
                                    eldispiter(index) = displ_im1((conn(i,j)-1)*numDOFPerNode + k);
                                    if(strcmp(iterationType,'NR'))
                                          eldisp(index) = displ_im1((conn(i,j)-1)*numDOFPerNode + k);
                                          eldispdot(index) = displdot_im1((conn(i,j)-1)*numDOFPerNode + k);
                                          eldispddot(index) = displddot_im1((conn(i,j)-1)*numDOFPerNode + k);
                                    elseif(strcmp(iterationType,'DI')||strcmp(iterationType,'LINEAR'))
                                          eldisp(index) = disp_s((conn(i,j)-1)*numDOFPerNode + k);
                                          eldispdot(index) = dispdot_s((conn(i,j)-1)*numDOFPerNode + k);
                                          eldispddot(index) = dispddot_s((conn(i,j)-1)*numDOFPerNode + k);
                                    end
                              end
                              index = index + 1;
                        end
                  end

                  %get concentrated terms associated with elemetn
                  [massConc,stiffConc,loadConc,model.joint,nodalTerms.concMass,nodalTerms.concStiff] = ConcMassAssociatedWithElement(conn(i,:),model.joint,nodalTerms.concMass,nodalTerms.concStiff,nodalTerms.concLoad);

                  elInput.concMass = massConc;
                  elInput.concStiff = stiffConc;
                  elInput.concLoad = loadConc;
                  elInput.disp = eldisp;

                  % specific to 'TD', but must be declared
                  elInput.dispm1= eldisp_sm1;

                  % specific to 'TNB' , but must be declared
                  elInput.dispdot = eldispdot;
                  elInput.dispddot = eldispddot;

                  elInput.x = elx;
                  elInput.y = ely;
                  elInput.z = elz;

                  if(el.rotationalEffects)
                        elInput.Omega = Omega;
                        elInput.OmegaDot = OmegaDot;
                  else
                        elInput.Omega = 0.0;
                        elInput.OmegaDot = 0.0;
                  end

                  elInput.displ_iter = eldispiter;

                  [elOutput] = calculateTimoshenkoElementNL(elInput,elStorage(i)); %calculate timoshenko element

                  conin = conn(i,:);


                  [Kg,Fg] = assembly(elOutput.Ke,elOutput.Fe,conn(i,:),numNodesPerEl,numDOFPerNode,Kg,Fg); %assemble element stiffness matrix and force vector

                  %         Erestotal = Erestotal + elOutput.Eres;
                  %................................................
            end %for
            %------- end element calculation and assembly ------------------

            %%
            %----------------------------------------------------------------------

            %%
            %Apply external loads to structure
            for i=1:length(Fexternal)
                  if(strcmp(analysisType,'TD'))
                        Fg(Fdof(i)) = Fg(Fdof(i)) + Fexternal(i)*delta_t^2;
                  end
                  if(strcmp(analysisType,'TNB'))
                        Fg(Fdof(i)) = Fg(Fdof(i)) + Fexternal(i);
                  end
            end

            %------ apply constraints on system -----------------------------------
            [Kg] = applyConstraints(Kg,model.jointTransform);
            [Fg] = applyConstraintsVec(Fg,model.jointTransform);





            %----------------------------------------------------------------------
            %%

            %Apply BCs to global system
            [KgTotal,FgTotal] = applyBC(Kg,Fg,BC,numDOFPerNode);

            solution = KgTotal\FgTotal;  %solve for displacements

            solution = model.jointTransform*solution; %transform to full dof listing

            if(model.nlOn)  %calculate norm between current iteration and last iteration
                  if(strcmp(iterationType,'NR'))
                        [unorm] = calcUnorm(displ_im1+solution,displ_im1);
                  else
                        [unorm] = calcUnorm(solution,displ_im1);
                  end
            else
                  unorm = 0.0;
            end

            if(strcmp(iterationType,'NR'))
                  %if newton raphson update u, udot, uddot at each iteration
                  displ_im1 = displ_im1 + solution;
                  cap_delta_displ = displ_im1 - dispData.displ_s;
                  displddot_im1 = timeInt.a3*(cap_delta_displ) - timeInt.a4*dispData.displdot_s - timeInt.a5*dispData.displddot_s;
                  displdot_im1  = -timeInt.a7*dispData.displdot_s -timeInt.a8*dispData.displddot_s + timeInt.a6*(cap_delta_displ);
            elseif(strcmp(iterationType,'DI')||strcmp(iterationType,'LINEAR'))
                  displ_im1 = solution;
            else
                  error('iteration type not supported, choose another')
            end

            iterationCount = iterationCount + 1;
      end %While

      %Calculate reaction at turbine base (hardwired to node number 1)
      reactionNodeNumber = platformTurbineConnectionNodeNumber;


      [FReaction] = calculateReactionForceAtNode(reactionNodeNumber,model,...
      mesh,el,elStorage,timeInt,dispData,displ_im1,rbData,Omega,OmegaDot,CN2H);

      %Calculate strain
      [elStrain] = calculateStrainForElements(numEl,numNodesPerEl,numDOFPerNode,conn,elementOrder,el,displ_im1,model.nlOn);

      % Size of eps and gam are mesh.numEl x 4
      eps_xx_0 = zeros(1,300);
      eps_xx_z = zeros(1,300);
      eps_xx_y = zeros(1,300);
      gam_xz_0 = zeros(1,300);
      gam_xz_y = zeros(1,300);
      gam_xy_0 = zeros(1,300);
      gam_xy_z = zeros(1,300);

      for jj = 1:length(elStrain)
            for ii = 1:4
                  eps_xx_0(jj*4+ii-4) = elStrain(jj).eps_xx_0(ii);
                  eps_xx_z(jj*4+ii-4) = elStrain(jj).eps_xx_z(ii);
                  eps_xx_y(jj*4+ii-4) = elStrain(jj).eps_xx_y(ii);
                  gam_xz_0(jj*4+ii-4) = elStrain(jj).gam_xz_0(ii);
                  gam_xz_y(jj*4+ii-4) = elStrain(jj).gam_xz_y(ii);
                  gam_xy_0(jj*4+ii-4) = elStrain(jj).gam_xy_0(ii);
                  gam_xy_z(jj*4+ii-4) = elStrain(jj).gam_xy_z(ii);
            end
      end

      if(iterationCount>=maxIterations)
            error('Maximum iterations exceeded.');
      end

      FReaction_sp1 =FReaction;
      displ_sp1 = displ_im1;
      % dispOut.displ_sp1 = displ_sp1;  %store displacement vector in dispOut

      % Specific to TNB, but must be declared
      displddot_sp1 = timeInt.a3*(displ_sp1-dispData.displ_s) - timeInt.a4*dispData.displdot_s - timeInt.a5*dispData.displddot_s; %store velocity vector in dispOut
      % dispOut.displddot_sp1 = displddot_sp1;
      displdot_sp1 = dispData.displdot_s + timeInt.a2*dispData.displddot_s + timeInt.a1*displddot_sp1;                    %store acceleration vector in dispOut

      % fprintf('%1.1e\n', FReaction_sp1(1));

end


function [Kg] = applyConstraints(Kg,transMatrix)
      %This function transforms a matrix by the transformation matrix to
      %enforce joint constraints
      % Kg = transMatrix'*(Kg*transMatrix);
      Kg = full(sparse(transMatrix')*sparse(Kg)*sparse(transMatrix));
end

function [Fg] = applyConstraintsVec(Fg,transMatrix)
      %This function transforms a vector by the transformation matrix to
      %enforce joint constraints
      Fg = transMatrix'*Fg;
end

function [unorm] = calcUnorm(unew,uold)
      %This function calculates a relative norm between two vectors: unew and
      %uold
      unorm = norm(unew-uold)/norm(unew);
end
