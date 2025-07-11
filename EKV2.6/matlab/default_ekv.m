function model=default_ekv();
model=struct('cox',3.45e-3,'xj',.15e-6,'vto',.7,'gamma',.7,'phi',.5,'kp',150e-6,...
    'e0',200e6,'ucrit',2.3e6,'dl',-7.5e-8,'dw',-7.5e-8,'lambda',0.8,'leta',0.3,'weta',0.2,...
    'q0',230e-6,'lk',.4e-6,'iba',2.0e8,'ibb',2e8,'ibn',.6,'theta',0,'ekvint',0,...
    'tcv',1e-3,'bex',-1.5,'ucex',.8,'ibbt',9e-4,'avto',0,'akp',0,'agamma',0,...
    'kf',0,'af',0,'nqs',0,'satlim',exp(4),'xqc',0.4,'ekvdyn',0,'type',0,'update',26,...
    'rsh',0,'hdif',0);