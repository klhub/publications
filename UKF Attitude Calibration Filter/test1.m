qa1 = [2 3 5 1]';
qa1 = qa1/norm(qa1);
qa2 = [1 2 3 4]'; 
qa2 = qa2/norm(qa2);
qa3 = [2 3 3 6]';
qa3 = qa3/norm(qa3);
qa  = [qa1 qa2 qa3]

qb1 = [-2 -3 -5 1]';
qb1 = qb1/norm(qb1);
qb2 = [-1 -2 -3 4]'; 
qb2 = qb2/norm(qb2);
qb3 = [-2 -3 -3 6]';
qb3 = qb3/norm(qb3);
qb  = [qb1 qb2 qb3]

%a=q_mult(qa(:,1),qb(:,1))

b=q_mult(qa,qb)
    