% test a "tower" with a mass ontop on a cart
% This test demonstrates that
% * the apparent stiffness of a beam is reduced by its own weight in longitudinal direction
% * the apparent stiffness of a beam is reduced by a mass ontop
% * the effective inertia of a beam is influenced by the moment of inertia of an attached mass

elastic_body_system= MultiBodySystem('elastic_body_system');
elastic_body_system.addGeneralizedCoordinate('x_cart')
elastic_body_system.addGeneralizedCoordinate('q_tow')

elastic_body_system.addInput('F_cart')
elastic_body_system.addParameter('g', [], 9.81);
elastic_body_system.addParameter('m_cart', [], 1000);
elastic_body_system.addParameter('m_top', [], 1000);
elastic_body_system.addParameter('I_top', [], 1000);
elastic_body_system.gravity(3) = -elastic_body_system.params.g;

cart= RigidBody('cart', [], elastic_body_system.params.m_cart);
cart.translate([elastic_body_system.dof.x_cart 0 0]);
elastic_body_system.addChild(cart)

load('tow_sid.mat')
tw_sid_ = SID(tw_sid, -1e-6, 'tow', elastic_body_system);
tower= ElasticBody(elastic_body_system.dof.q_tow, tw_sid_, 'tower');
cart.addChild(tower)

top_mass= RigidBody('top_mass', [], elastic_body_system.params.m_top);
top_mass.I(2, 2)= elastic_body_system.params.I_top;
tower.addChild(top_mass)


cart.applyForce([elastic_body_system.inputs.F_cart; 0; 0])
elastic_body_system.addOutput('x', elastic_body_system.dof.x_cart)
elastic_body_system.addOutput('x_tow', elastic_body_system.dof.q_tow)
elastic_body_system.completeSetup;

eom = elastic_body_system.getEOM


%%
elastic_body_system.eomDae('testElasticBody_DAE.m')

t_end = 50;
params = elastic_body_system.getParameterArray();
f_min = 0.01;
f_max = 5;
excitation = @(t)chirp(t, f_min, t_end, f_max);
F = @(t,Y,YP) testElasticBody_DAE(t, Y, YP, [excitation(t)], [], params);

x0 = [0 0 0 0]';
dx0= [x0(elastic_body_system.getNumDOF+1:end)' 0*x0(elastic_body_system.getNumDOF+1:end)']';

opt= odeset(AbsTol=1e-5, RelTol=1e-3, MaxStep=0.01);

sim_res= ode15i(F, [0,t_end], x0, dx0, opt);

ff = interp1([sim_res.x(1) sim_res.x(end)], [f_min f_max], sim_res.x);
plot(ff, sim_res.y(elastic_body_system.dof_idx.q_tow, :))
grid on