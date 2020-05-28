%piano.m
clear all

%define rule, EDIT THIS
zero = 0; %000
one = 1;  %001
two = 1;  %010
three = 1;%011
four = 0; %100
five = 1; %101
six = 1;  %110
seven = 0;%111

%starting generation, EDIT THIS
current = [1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0];

%number of generations, EDIT THIS
numgen = 20;

%-----------------------------------------------------------------------

rule = [seven, six, five, four, three, two, one, zero];
ruledata = ["111", "110", "101", "100", "011", "010", "001", "000"];

%specify the music to be played:
nmax= 0; %number of notes
currgen = 0;
tnote = []; %array of onset times of the notes (s)
inote = []; %array of string indices of the notes

bw = [];

for l=1:numgen %number of generations
    new = zeros(1, 13);
    
    for i=1:13
        curr = i;
        prev = i-1;
        next = i+1;
        
        if (prev == 0)
            prev = 13;
        end
        if (next == 14)
            next = 1;
        end
        
        test = string(current(prev)) + string(current(curr)) + string(current(next));
        
        for r=1:8
            if (strcmp(test, ruledata(r)))
                new(curr) = rule(r);
                break;
            end
        end
    end
    
    %add to arrays
    for g=1:13
        if(new(g) == 1)
            nmax = nmax + 1;
            tnote(length(tnote)+1) = currgen;
            inote(length(inote)+1) = g;
        end
    end
    
    currgen = currgen + 1;
    
    bw(currgen+1, 1:13) = current;
    
    current = new;
end

dnote = ones(1, length(tnote));
anote = ones(1, length(tnote));

%----------------------------------------------------------------------

%initialize string parameters:
L=1;          % length of all strings
J=81;dx=L/(J-1); %number of points/string, space step
flow = 220; %frequency of string with lowest pitch (1/s)
nstrings = 13; %number of strings
for i=1:nstrings
  f(i)=flow*2^((i-1)/12);  % frequency (1/s)
  tau(i)=1.2*(440/f(i));   % decay time (s)
  M(i)=1;                  % mass/length
  T(i)=M(i)*(2*L*f(i))^2;  % tension
  R(i)=(2*M(i)*L^2)/(tau(i)*pi^2); % damping constant 
  %Find the largest stable timestep for string i:
  dtmax(i) = - R(i)/T(i) + sqrt((R(i)/T(i))^2 + dx^2/(T(i)/M(i)));
end
%The timestep of the computation has to be stable for all strings:
dtmaxmin = min(dtmax);
%Now set dt and nskip such that:
%dt<=dtmaxmin, nskip is a positive integer, and dt*nskip = 1/8192.
%Also, make nskip as small as possible, given the above criteria.
nskip = ceil(1/(8192*dtmaxmin));
dt=1/(8192*nskip);
tmax=tnote(nmax)+dnote(nmax); clockmax=ceil(tmax/dt);
%initialize an array that will be used to tell a string 
%when to stop vibrating:
tstop=zeros(nstrings,1);
%initialize arrays to store state of the strings:
H=zeros(nstrings,J);
V=zeros(nstrings,J);

%(xh1,xh2)=part of string hit by hammer:
xh1=0.25*L;xh2=0.35*L;
%list of points hit by hammer:
jstrike=ceil(1+xh1/dx):floor(1+xh2/dx);
j=2:(J-1); %list of interior points
%initialize array to store soundwave:
count=0; %initialize count of samples of recorded soundwave
S=zeros(1,ceil(clockmax/nskip)); %array to record soundwave
tsave = zeros(1,ceil(clockmax/nskip)); %array for times of samples

n=1 ; %initialize note counter
for clock=1:clockmax
  t=clock*dt;
  while((n<=nmax) && tnote(n)<=t)
    V(inote(n),jstrike)=anote(n); %strike string inote(n)
                                  %with amplitude anote(n)
    tstop(inote(n))=t+dnote(n);   %record future stop time
    n=n+1;                        %increment note counter
  end
  for i=1:nstrings
    if(t > tstop(i))
      H(i,:)=zeros(1,J);
      V(i,:)=zeros(1,J);
    else
      V(i,j)=V(i,j) ... %solving 2 eqs on pg 6
            +(dt/dx^2)*(T(i)/M(i))*(H(i,j+1)-2*H(i,j)+H(i,j-1)) ...
            +(dt/dx^2)*(R(i)/M(i))*(V(i,j+1)-2*V(i,j)+V(i,j-1));
      H(i,j)=H(i,j)+dt*V(i,j); %H(x, t + dealta t)
    end
  end
  if(mod(clock,nskip)==0)
    count=count+1;
S(count)=sum(H(:,2)); %sample the sound at the present time
    tsave(count)=t;      %record the time of the sample
  end
end
soundsc(S(1:count)) %listen to the sound
imshow(bw, 'InitialMagnification', 'fit');


