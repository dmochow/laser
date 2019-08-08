myArray=zeros(10,10); % you need to figure out the size or your grid
values=randn(100,1); % you need to figure 

currentPointer=1; % where you are in the values
for row=1:10
    myArray(row,:)=values(currentPointer:currentPointer+10-1);
    currentPointer=currentPointer+10;
end
