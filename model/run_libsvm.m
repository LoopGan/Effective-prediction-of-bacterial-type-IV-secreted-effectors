train_data=xlsread('train_data.xlsx');                                                                                   
train_label=xlsread('train_label.xlsx');
test_data=xlsread('test_data.xlsx');
test_label=xlsread('test_label.xlsx');
mapminmax1                                                                    
[bestacc,bestc,bestg] = SVMcgForClass(train_label,train_data,-8,8,-8,8,5,1,1,0.5); 
cmd = ['-c ', num2str(bestc), ' -g ', num2str(bestg)];
model = svmtrain(train_label, train_data, cmd);
acc=svmpredict(test_label,test_data,model); 
[predicted_label,accuracy,score]=svmpredict(test_label,test_data,model); 
[SE,SP,ACC,MCC]=Runget_FTP_MCC(test_label,predicted_label)

