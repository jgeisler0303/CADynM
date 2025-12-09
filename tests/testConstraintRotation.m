% function testConstraintRotation
% Simple test to verify constraint forces
testSystem= MultiBodySystem('rot_constr_test');
testSystem.addGeneralizedCoordinate('phi')
testSystem.addConstraintCoordinate('Fc_x')
testSystem.addConstraintCoordinate('Fc_y')
testSystem.addInput('Fx')
testSystem.addParameter('m');
testSystem.addParameter('l');
testSystem.addParameter('g');
testSystem.gravity(2) = -testSystem.params.g;

rotatingBody= RigidBody('rotor', 'Rotating body with off-center mass');
rotatingBody.m= testSystem.params.m;

rotatingBody.translate([testSystem.doc.Fc_x testSystem.doc.Fc_y 0]);
rotatingBody.rotateLocalAxis('z', testSystem.dof.phi);
rotatingBody.translate([testSystem.params.l 0 0]);

testSystem.addChild(rotatingBody)

testSystem.applyForce([testSystem.inputs.Fx; 0; 0], rotatingBody)

testSystem.prepareKinematics()
% testSystem.eomDae('testConstraintRotation_DAE.m')
Fz = testSystem.getConstraintForces()