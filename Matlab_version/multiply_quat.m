function q = multiply_quat(q1,q2)
    q = zeros(4,1);

    qv1 = q1(1:3);
    qv2 = q2(1:3);
    qs1 = q1(4);
    qs2 = q2(4);
    q(1:3) = qv1*qs2 + qv2*qs1 - cross(qv1,qv2);
    q(4) = qs1*qs2 - dot(qv1,qv2);
end
