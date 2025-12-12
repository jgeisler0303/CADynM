% function testPendulum
% The classic pendulum on cart model
pendulum_model= MultiBodySystem('pendulum');
pendulum_model.addGeneralizedCoordinate('x_cart')
pendulum_model.addGeneralizedCoordinate('phi_pend')

pendulum_model.addInput('F_cart')
pendulum_model.addParameter('g');
pendulum_model.addParameter('m_cart');
pendulum_model.addParameter('m_pend');
pendulum_model.addParameter('l');
pendulum_model.gravity(3) = -pendulum_model.params.g;

cart= RigidBody('cart');
cart.m= pendulum_model.params.m_cart;
cart.translate([pendulum_model.dof.x_cart 0 0]);
pendulum_model.addChild(cart)

pendulum= RigidBody('pendulum');
pendulum.m= pendulum_model.params.m_pend;
pendulum.rotateLocalAxis('y', pendulum_model.dof.phi_pend);
pendulum.translate([0 0 pendulum_model.params.l]);
cart.addChild(pendulum)

pendulum_model.applyForce([pendulum_model.inputs.F_cart; 0; 0], cart)
pendulum_model.addOutput('x', pendulum_model.dof.x_cart)
pendulum_model.completeSetup;

[~,~]=mkdir('generated');
matlabTemplateEngine('generated/model_parameters.m', 'model_parameters.m.mte', pendulum_model)
matlabTemplateEngine('generated/model_indices.m', 'model_indices.m.mte', pendulum_model)
matlabTemplateEngine('generated/param.hpp', 'param.hpp.mte', pendulum_model)
matlabTemplateEngine('generated/direct.hpp', 'direct.hpp.mte', pendulum_model)
% pendulum_model.prepareKinematics()
% pendulum_model.eomDae('testPendulum.m')
% eom = pendulum_model.getEOM()