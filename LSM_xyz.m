function [UserPos] = LSM_xyz(SatPos, D, InitPos)
MaxErr = 1e-3;
MaxIter = 20;
stopflag = 1;
k=0;
UserPos = InitPos;
N = length(InitPos);
rho = zeros(length(D),1);
H = zeros(length(D),N);
Niter = 0;
while stopflag
    k=k+1;
    Niter = Niter + 1;
    for m = 1:length(D)
        rho(m) = norm(SatPos(1:N,m)-UserPos(1:N));
        for n = 1:N
            H(m,n) = -(SatPos(n,m)-UserPos(n))/norm(SatPos(1:N,m)-UserPos(1:N));
        end
    end
    DU=UserPos; %���������� ��� ������� ��������
    dPR=D-rho; %������� �� ����������
    UserPos=UserPos+((H'*H)^(-1))*H'*dPR; %��������� �������� �������
    dX=UserPos-DU; %���������� ������
    DOP=sqrt(trace((H'*H)^(-1))); %�������������� ������
    %�������� ��������
    if norm(dX)<=MaxErr
        stopflag=0; %����� �� �����, ���� ���������� ��������� ��������
    elseif Niter>=MaxIter
            stopflag=0; %����� �� �����, ���� ���������� ���������� ����� ��������
    end
end
end