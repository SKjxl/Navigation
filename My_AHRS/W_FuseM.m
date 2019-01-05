function [...
    nextStates, ... % state vector after fusion of measurements
    nextP...
    Tbn,zk] ... % state covariance matrix after fusion of corrections
    = W_FuseM( ...
    states, ... % predicted state
    P, ... % predicted covariance
    Zk, ...
    Tbn,...
    initMag)  % Estimated coordinate transformation matrix from body to NED frame

q0 = states(1);
q1 = states(2);
q2 = states(3);
q3 = states(4);

mx = initMag(1);
my = initMag(2);
mz = initMag(3);
R = 9e-6*eye(3);
H_mag = W_calcHm(mx,my,mz,q0,q1,q2,q3);
varInnov = (H_mag*P*H_mag' + R);
Kfusion = (P*H_mag')*inv(varInnov);

% Calculate the predicted magnetic declination
% z = [Tbn'*[0;0;1];Tbn'*[mx;my;mz]];
zk = Tbn'*[mx my mz]';
% Calculate the measurement innovation
innovation = Zk - zk ;

% correct the state vector
% states(1:3) = 0;
states = states+ Kfusion * innovation;

% the first 3 states represent the angular misalignment vector. This is
% is used to correct the estimate quaternion
% Convert the error rotation vector to its equivalent quaternion
% error = truth - estimate
% rotationMag = sqrt(states(1)^2 + states(2)^2 + states(3)^2);
% if rotationMag<1e-6
%     deltaQuat = single([1;0;0;0]);
% else
%     deltaQuat = [cos(0.5*rotationMag); [states(1);states(2);states(3)]/rotationMag*sin(0.5*rotationMag)];
% end
% 
% % Update the quaternion states by rotating from the previous attitude through
% % the delta angle rotation quaternion
% nextQuat = [quat(1)*deltaQuat(1)-transpose(quat(2:4))*deltaQuat(2:4); quat(1)*deltaQuat(2:4) + deltaQuat(1)*quat(2:4) + cross(quat(2:4),deltaQuat(2:4))];
quat = states(1:4);
Tbn = Quat2Tbn(quat);
% normalise the updated quaternion states
quatMag = sqrt(quat(1)^2 + quat(2)^2 + quat(3)^2 + quat(4)^2);
if (quatMag > 1e-6)
    quat = quat / quatMag;
end

states(1:4) = quat;

% correct the covariance P = P - K*H*P
P = P - Kfusion*H_mag*P;

% Force symmetry on the covariance matrix to prevent ill-conditioning
% of the matrix which would cause the filter to blow-up
P = 0.5*(P + transpose(P));

% % ensure diagonals are positive
% for i=1:9
%     if P(i,i) < 0
%         P(i,i) = 0;
%     end
% end

% Set default output for states and covariance
nextP = P;
nextStates = states;

end