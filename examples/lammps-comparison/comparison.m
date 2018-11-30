qs=ncread('../MD0.nc','qs');
qs2lmp=qs*0.01865707;

for i=1:1000
    lmpq(:,i)=reshape(transpose(importfile(['dump',num2str(i-1),'.xyz'])),18,1)-reshape(transpose(importfile('dump0.xyz')),18,1)
end

figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'on');
plot1 = plot(1:1000,lmpq(15,:),1:1000,qs2lmp(15,:),".");
set(plot1(1),'DisplayName','data1 from LAMMPS');
set(plot1(2),'DisplayName','data2 from MD','Marker','.','LineStyle','none');
ylabel({'\Delta q(Angstroms)'});
xlabel({'MD times'});
box(axes1,'on');
legend(axes1,'show');