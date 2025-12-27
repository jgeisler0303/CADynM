% function testExternalLinear
% Simple test to verify constraint forces
testSystem= MultiBodySystem('lin_external_test', 'x', 'u');
testSystem.addExternal('F', [testSystem.dof.x, testSystem.inputs.u, testSystem.getTimeDeriv(testSystem.dof.x)])
testSystem.addParameter({'m' 'l'});

movingBody= RigidBody('mover', 'Simple mass moving');
movingBody.m= testSystem.params.m;

movingBody.translate([testSystem.dof.x+testSystem.params.l 0 0]);

testSystem.addChild(movingBody)

movingBody.applyForce([testSystem.externals.F; 0; 0])

testSystem.completeSetup()
[~,~]=mkdir('generated');
testSystem.eomDae('generated/testExternalLinear_DAE.m')
eom = testSystem.getEOM()
testSystem.getExternalDerivs(1)

testSystem.replaceVars(testSystem.getK)
testSystem.replaceVars(testSystem.getC)
testSystem.replaceVars(testSystem.getB)