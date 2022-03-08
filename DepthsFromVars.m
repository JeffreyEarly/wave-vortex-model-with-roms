function z = DepthsFromVars(Vtransform, hc, sc_r, Cs_r, sc_w, Cs_w, h, igrid, idims, zeta)
%% DepthsFromVars
%
% March 7th, 2022 --- Jeffrey J. Early
%
% Taken from roms-matlab/utility/depths.m, svn $Id: depths.m 1102 2022-01-06 21:50:44Z arango $
%
% history file, e.g. niskine_his_00040.nc
% Vtransform:long_name = "vertical terrain-following transformation equation" ;
% hc:long_name = "S-coordinate parameter, critical depth"
% s_rho:long_name = "S-coordinate at RHO-points" ;
% Cs_r:long_name = "S-coordinate stretching curves at RHO-points" ;
% s_w:long_name = "S-coordinate at W-points" ;
% Cs_w:long_name = "S-coordinate stretching curves at W-points" ;
%
% geometry file, e.g., NISKINEthilo_4km_negLons.nc
% h:long_name = "Final bathymetry at RHO-points" ;
%
% history2 file, e.g. niskine_his2_00040.nc
% zeta


N=length(sc_r);
Np=N+1;

if (length(sc_w) == N),
    sc_w=[-1 sc_w'];
    Cs_w=[-1 Cs_w'];
end

%------------------------------------------------------------------------
% Get bottom topography.
%------------------------------------------------------------------------

[Lp, Mp]=size(h);
L=Lp-1;
M=Mp-1;

switch ( igrid )
    case 1
        if (idims), h=h'; end
    case 2
        hp=0.25.*(h(1:L,1:M)+h(2:Lp,1:M)+h(1:L,2:Mp)+h(2:Lp,2:Mp));
        if (idims), hp=hp'; end
    case 3
        hu=0.5.*(h(1:L,1:Mp)+h(2:Lp,1:Mp));
        if (idims), hu=hu'; end
    case 4
        hv=0.5.*(h(1:Lp,1:M)+h(1:Lp,2:Mp));
        if (idims), hv=hv'; end
    case 5
        if (idims), h=h'; end
end

% Set critical depth parameter.
if isempty(hc)
    hc=min(min(h));
end

%------------------------------------------------------------------------
% Get free-surface
%------------------------------------------------------------------------

if isempty(zeta)
    zeta=zeros([Lp Mp]);
end

switch ( igrid )
    case 1
        if (idims), zeta=zeta'; end
    case 2
        zetap=0.25.*(zeta(1:L,1:M )+zeta(2:Lp,1:M )+                      ...
            zeta(1:L,2:Mp)+zeta(2:Lp,2:Mp));
        if (idims), zetap=zetap'; end
    case 3
        zetau=0.5.*(zeta(1:L,1:Mp)+zeta(2:Lp,1:Mp));
        if (idims), zetau=zetau'; end
    case 4
        zetav=0.5.*(zeta(1:Lp,1:M)+zeta(1:Lp,2:Mp));
        if (idims), zetav=zetav'; end
    case 5
        if (idims), zeta=zeta'; end
end

%------------------------------------------------------------------------
% Compute depths.
%------------------------------------------------------------------------

if (Vtransform == 1)
    switch ( igrid )
        case 1
            for k=1:N
                z0=(sc_r(k)-Cs_r(k))*hc + Cs_r(k).*h;
                z(:,:,k)=z0 + zeta.*(1.0 + z0./h);
            end
        case 2
            for k=1:N
                z0=(sc_r(k)-Cs_r(k))*hc + Cs_r(k).*hp;
                z(:,:,k)=z0 + zetap.*(1.0 + z0./hp);
            end
        case 3
            for k=1:N
                z0=(sc_r(k)-Cs_r(k))*hc + Cs_r(k).*hu;
                z(:,:,k)=z0 + zetau.*(1.0 + z0./hu);
            end
        case 4
            for k=1:N
                z0=(sc_r(k)-Cs_r(k))*hc + Cs_r(k).*hv;
                z(:,:,k)=z0 + zetav.*(1.0 + z0./hv);
            end
        case 5
            z(:,:,1)=-h;
            for k=2:Np
                z0=(sc_w(k)-Cs_w(k))*hc + Cs_w(k).*h;
                z(:,:,k)=z0 + zeta.*(1.0 + z0./h);
            end
    end
elseif (Vtransform == 2)
    switch ( igrid )
        case 1
            for k=1:N
                z0=(hc.*sc_r(k)+Cs_r(k).*h)./(hc+h);
                z(:,:,k)=zeta+(zeta+h).*z0;
            end
        case 2
            for k=1:N
                z0=(hc.*sc_r(k)+Cs_r(k).*hp)./(hc+hp);
                z(:,:,k)=zetap+(zetap+hp).*z0;
            end
        case 3
            for k=1:N
                z0=(hc.*sc_r(k)+Cs_r(k).*hu)./(hc+hu);
                z(:,:,k)=zetau+(zetau+hu).*z0;
            end
        case 4
            for k=1:N
                z0=(hc.*sc_r(k)+Cs_r(k).*hv)./(hc+hv);
                z(:,:,k)=zetav+(zetav+hv).*z0;
            end
        case 5
            for k=1:Np
                z0=(hc.*sc_w(k)+Cs_w(k).*h)./(hc+h);
                z(:,:,k)=zeta+(zeta+h).*z0;
            end
    end
end

end