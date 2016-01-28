function Fm=FlocofFmac(Fmac,hoR,alpha)

poly=[hoR*Fmac*alpha -1 (2*Fmac+hoR*Fmac)];

rts=roots(poly);

Fm=min(rts);
end
