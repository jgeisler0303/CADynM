% test a "tower" with a mass ontop on a cart
% This test demonstrates that
% * the apparent stiffness of a beam is reduced by its own weight in longitudinal direction
% * the apparent stiffness of a beam is reduced by a mass ontop
% * the effective inertia of a beam is influenced by the moment of inertia of an attached mass

elastic_body_system= MultiBodySystem('elastic_body_system');
elastic_body_system.addGeneralizedCoordinate('x_cart')
elastic_body_system.addGeneralizedCoordinate('q_tow')

elastic_body_system.addInput('F_cart')
elastic_body_system.addParameter('g');
elastic_body_system.addParameter('m_cart');
elastic_body_system.addParameter('m_top');
elastic_body_system.addParameter('I_top');
elastic_body_system.gravity(3) = -elastic_body_system.params.g;

cart= RigidBody('cart');
cart.m= elastic_body_system.params.m_cart;
cart.translate([elastic_body_system.dof.x_cart 0 0]);
elastic_body_system.addChild(cart)

load('tow_sid.mat')
tw_sid_ = SID(tw_sid, -1e-6, 'tow', elastic_body_system);
tower= ElasticBody(elastic_body_system.dof.q_tow, tw_sid_, 'tower');
cart.addChild(tower)

top_mass= RigidBody('top_mass');
top_mass.m= elastic_body_system.params.m_top;
top_mass.I(2, 2)= elastic_body_system.params.I_top;
tower.addChild(top_mass)


elastic_body_system.applyForce([elastic_body_system.inputs.F_cart; 0; 0], cart)
elastic_body_system.addOutput('x', elastic_body_system.dof.x_cart)
elastic_body_system.addOutput('x_tow', elastic_body_system.dof.q_tow)
elastic_body_system.completeSetup;

eom = elastic_body_system.getEOM
