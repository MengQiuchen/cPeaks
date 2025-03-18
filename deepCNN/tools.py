import torch
import numpy as np
import torch.nn as nn
from tqdm import tqdm
import matplotlib.pyplot as plt
from torch.utils.data import TensorDataset
from torch.utils.data import Dataset, DataLoader
from torch.utils.data import random_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score


def seq2tensor(seq):
    num_seq = [0 if base == 'A' else 1 if base == 'C' else 2 if base == 'G' else 3 if base == 'T' else 4 for base in seq]

    # 进行onehot编码
    onehot_seq = torch.zeros((5, len(num_seq)))
    for i, num in enumerate(num_seq):
        onehot_seq[num, i] = 1
    
    onehot_seq = onehot_seq[:4,:]
    
    # 统一序列长度为2000，并进行填充
    if len(num_seq) >= 2000:
        padded_seq = onehot_seq[:, :2000]
    else:
        padded_seq = torch.zeros((4, 2000))
        padded_seq[:, :len(num_seq)] = onehot_seq
        
    return padded_seq
        
class ChromosomeDataset_Train(Dataset):
    
    def __init__(self, pos_fa_file,neg_fa_file):

        # 从fa文件读取数据
        with open(pos_fa_file, 'r') as file:
            tmp = [line.strip() for  line in file]
            bed1  = [line for i, line in enumerate(tmp) if i % 2 == 0]
            data1 = [line.upper() for i, line in enumerate(tmp) if i % 2 == 1]
            
        with open(neg_fa_file, 'r') as file:
            tmp = [line.strip() for  line in file]
            bed2  = [line for i, line in enumerate(tmp) if i % 2 == 0]
            data2 = [line.upper() for i, line in enumerate(tmp) if i % 2 == 1]
            
        self.data = data1+data2
        self.bed = bed1+bed2
        
        self.label = [1]*len(data1)+[0]*len(data2)

        print('fasta loaded!')

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sequence = self.data[idx]
        sequence_tensor = seq2tensor(sequence)
        return self.bed[idx],sequence_tensor,self.label[idx]

class ChromosomeDataset_Test(Dataset):
    
    def __init__(self, pos_fa_file):

        # 从fa文件读取数据
        with open(pos_fa_file, 'r') as file:
            tmp = [line.strip() for  line in file]
            bed1  = [line for i, line in enumerate(tmp) if i % 2 == 0]
            data1 = [line.upper() for i, line in enumerate(tmp) if i % 2 == 1]
            
       
        self.data = data1
        self.bed = bed1
        
        print('fasta loaded!')

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sequence = self.data[idx]
        sequence_tensor = seq2tensor(sequence)
        return self.bed[idx],sequence_tensor
        

def train(model, train_loader, optimizer, criterion, device):
    model.train()  # 将模型设置为训练模式
    model.to(device)
    
    total_loss = 0.0  # 总损失
    outputs_list = []
    targets_list = []
    for batch_idx, (_, inputs, targets) in enumerate(train_loader):

        inputs, targets = inputs.to(device), targets.to(device)

        optimizer.zero_grad()  # 梯度清零

        outputs = model(inputs)
        loss = criterion(outputs, targets)  # 计算损失
        loss.backward()  # 反向传播
        optimizer.step()  # 更新参数

        outputs_list.append(outputs)
        targets_list.append(targets)

        total_loss += loss.item()

    outputs_concat = torch.cat(outputs_list, dim=0)
    _, predicted_concat = outputs_concat.max(1)
    
    targets_concat = torch.cat(targets_list, dim=0).cpu().numpy()
    predicted_concat = predicted_concat.cpu().numpy()
    outputs_concat = outputs_concat[:,1].detach().cpu().numpy()

    avg_acc = accuracy_score(targets_concat, predicted_concat)  
    avg_precision = precision_score(targets_concat, predicted_concat, average='weighted')
    avg_recall = recall_score(targets_concat, predicted_concat, average='weighted')  
    avg_f1 = f1_score(targets_concat, predicted_concat, average='weighted')  
    avg_auc = roc_auc_score(targets_concat,outputs_concat, average='weighted')
    avg_loss = total_loss/len(train_loader)    

    print('Train Loss: {:.4f}, Acc: {:.4f}, Precision: {:.4f}, Recall: {:.4f}, F1: {:.4f}, AUC: {:.4f}'.format(
        avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc))
    
    return (avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc)

new_indices = torch.tensor([3, 2, 1, 0])     
def train_aug(model, train_loader, optimizer, criterion, device):

    model.train()  # 将模型设置为训练模式
    model.to(device)
    
    total_loss = 0.0  # 总损失
    outputs_list = []
    targets_list = []
    for _, inputs, targets in train_loader:

        inputs, targets = inputs.to(device), targets.unsqueeze(1).to(device)

        optimizer.zero_grad()  # 梯度清零

        outputs = model(inputs)
        
        one_hot = torch.zeros(targets.size(0), 2).to(device)
        one_hot.scatter_(1, targets, 1)
        loss = criterion(outputs, one_hot)
        
        loss.backward()  # 反向传播
        optimizer.step()  # 更新参数

        outputs_list.append(outputs)
        targets_list.append(targets)

        total_loss += loss.item()

        #数据增强训练
        optimizer.zero_grad()  # 梯度清零
        outputs = model(torch.flip(inputs, [2]).index_select(1, new_indices.to(device)))
        loss = criterion(outputs, one_hot)
        loss.backward()  # 反向传播
        optimizer.step()  # 更新参数
        
    predicted_concat = torch.cat(outputs_list, dim=0).detach().cpu().numpy()
    targets_concat = torch.cat(targets_list, dim=0).detach().cpu().numpy()

    avg_acc = accuracy_score(targets_concat, predicted_concat[:,1]>0.5)  
    avg_precision = precision_score(targets_concat, predicted_concat[:,1]>0.5, average='weighted')
    avg_recall = recall_score(targets_concat, predicted_concat[:,1]>0.5, average='weighted')  
    avg_f1 = f1_score(targets_concat, predicted_concat[:,1]>0.5, average='weighted')  
    avg_auc = roc_auc_score(targets_concat,predicted_concat[:,1], average='weighted')
    avg_loss = total_loss/len(train_loader)    

    print('Train Loss: {:.4f}, Acc: {:.4f}, Precision: {:.4f}, Recall: {:.4f}, F1: {:.4f}, AUC: {:.4f}'.format(
        avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc))
    
    return (avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc)
    
def test(model, test_loader, criterion, device):
    model.eval()  # 将模型设置为评估模式
    model.to(device)
    
    total_loss = 0.0  # 总损失
    outputs_list = []
    targets_list = []
    
    with torch.no_grad():  # 关闭梯度计算
        for _, inputs, targets in test_loader:
            
            inputs, targets = inputs.to(device), targets.unsqueeze(1).to(device)
            
            outputs = model(inputs)
            one_hot = torch.zeros(targets.size(0), 2).to(device)
            one_hot.scatter_(1, targets, 1)
            loss = criterion(outputs, one_hot)
            
            outputs_list.append(outputs)
            targets_list.append(targets)
            
            total_loss += loss.item()

    avg_loss = total_loss / len(test_loader)  # 平均损失
    
    outputs_concat = torch.cat(outputs_list, dim=0)
    _, predicted_concat = outputs_concat.max(1)
    targets_concat = torch.cat(targets_list, dim=0).cpu().numpy()
    predicted_concat = predicted_concat.cpu().numpy()
    outputs_concat = outputs_concat[:,1].detach().cpu().numpy()

    avg_acc = accuracy_score(targets_concat, predicted_concat)  
    avg_precision = precision_score(targets_concat, predicted_concat, average='weighted')
    avg_recall = recall_score(targets_concat, predicted_concat, average='weighted')  
    avg_f1 = f1_score(targets_concat, predicted_concat, average='weighted')  
    avg_auc = roc_auc_score(targets_concat,outputs_concat, average='weighted')

    
    print('Test Loss: {:.4f}, Acc: {:.4f}, Precision: {:.4f}, Recall: {:.4f}, F1: {:.4f}, AUC: {:.4f}'.format(
        avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc))
    
    return (avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc)

def plot_res(train_res,test_res,num_epochs,save_path=None):
    
    epoch_list = range(1, num_epochs+1)

    avg_loss_train = [res[0] for res in train_res]
    avg_acc_train = [res[1] for res in train_res]
    avg_precision_train = [res[2] for res in train_res]
    avg_recall_train = [res[3] for res in train_res]
    avg_f1_train = [res[4] for res in train_res]
    avg_auc_train = [res[5] for res in train_res]

    avg_loss_test = [res[0] for res in test_res]
    avg_acc_test = [res[1] for res in test_res]
    avg_precision_test = [res[2] for res in test_res]
    avg_recall_test = [res[3] for res in test_res]
    avg_f1_test = [res[4] for res in test_res]
    avg_auc_test = [res[5] for res in test_res]
    

    # 绘制图表
    plt.figure(figsize=(10, 6))

    # Average Loss
    plt.subplot(2, 3, 1)
    plt.plot(epoch_list, avg_loss_train, label='Training')
    plt.plot(epoch_list, avg_loss_test, label='Testing')
    plt.xlabel('Epoch')
    plt.ylabel('Average Loss')
    plt.title('Average Loss')
    plt.legend()

    # Average Accuracy
    plt.subplot(2, 3, 2)
    plt.plot(epoch_list, avg_acc_train, label='Training')
    plt.plot(epoch_list, avg_acc_test, label='Testing')
    plt.xlabel('Epoch')
    plt.ylabel('Average Accuracy')
    plt.title('Average Accuracy')
    plt.legend()

    # Average Precision
    plt.subplot(2, 3, 3)
    plt.plot(epoch_list, avg_precision_train, label='Training')
    plt.plot(epoch_list, avg_precision_test, label='Testing')
    plt.xlabel('Epoch')
    plt.ylabel('Average Precision')
    plt.title('Average Precision')
    plt.legend()

    # Average Recall
    plt.subplot(2, 3, 4)
    plt.plot(epoch_list, avg_recall_train, label='Training')
    plt.plot(epoch_list, avg_recall_test, label='Testing')
    plt.xlabel('Epoch')
    plt.ylabel('Average Recall')
    plt.title('Average Recall')
    plt.legend()

    # Average F1
    plt.subplot(2, 3, 5)
    plt.plot(epoch_list, avg_f1_train, label='Training')
    plt.plot(epoch_list, avg_f1_test, label='Testing')
    plt.xlabel('Epoch')
    plt.ylabel('Average F1')
    plt.title('Average F1')
    plt.legend()

    # Average AUC
    plt.subplot(2, 3, 6)
    plt.plot(epoch_list, avg_auc_train, label='Training')
    plt.plot(epoch_list, avg_auc_test, label='Testing')
    plt.xlabel('Epoch')
    plt.ylabel('Average AUC')
    plt.title('Average AUC')
    plt.legend()
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches='tight')
    plt.show()


class ChromosomeDataset_Train_shape(Dataset):
    def __init__(self, bed_file,fa_file):

        # 从fa文件读取数据
        with open(fa_file, 'r') as file:
            tmp = [line.strip() for  line in file]
            bed1  = [line for i, line in enumerate(tmp) if i % 2 == 0]
            data1 = [line.upper() for i, line in enumerate(tmp) if i % 2 == 1]
            
        with open(bed_file, 'r') as file:
            tmp = [line.strip().split('\t')[-2:] for  line in file]
            tmp = [[int(int(i[0])>170),int(int(i[1])>170)] for i in tmp]
    
            
        self.data = data1
        self.bed = bed1
        
        self.label = torch.tensor(tmp)

        print('fasta loaded!')

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sequence = self.data[idx]
        sequence_tensor = seq2tensor(sequence)
        return self.bed[idx],sequence_tensor,self.label[idx]
    
class ChromosomeDataset_Test_shape(Dataset):
    def __init__(self,fa_file):

        # 从fa文件读取数据
        with open(fa_file, 'r') as file:
            tmp = [line.strip() for  line in file]
            bed1  = [line for i, line in enumerate(tmp) if i % 2 == 0]
            data1 = [line.upper() for i, line in enumerate(tmp) if i % 2 == 1]
       
        self.data = data1
        self.bed = bed1
        
        print('fasta loaded!')

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        sequence = self.data[idx]
        sequence_tensor = seq2tensor(sequence)
        return self.bed[idx],sequence_tensor
    
def train_shape(model, train_loader, optimizer, criterion, device):
    model.train()
    total_loss = 0.0  # 总损失
    outputs_list = []
    targets_list = []
    for _, inputs, targets in train_loader:
        inputs, targets = inputs.to(device), targets.to(device)

        optimizer.zero_grad()  # 梯度清零

        outputs = model(inputs)
        loss = criterion(outputs, targets.float())  # 计算损失
        loss.backward()  # 反向传播
        optimizer.step()  # 更新参数

        outputs_list.append(outputs)
        targets_list.append(targets)

        total_loss += loss.item()
    
    avg_loss = total_loss / len(train_loader)
    
    return avg_loss
    

def test_shape(model, test_loader, optimizer, criterion, device):
    model = model.eval()
    total_loss = 0.0  # 总损失
    outputs_list = []
    targets_list = []
    
    with torch.no_grad():
        for _, inputs, targets in test_loader:
            inputs, targets = inputs.to(device), targets.to(device)

            outputs = model(inputs)
            loss = criterion(outputs, targets.float())  # 计算损失

            outputs_list.append(outputs)
            targets_list.append(targets)

            total_loss += loss.item()
        avg_loss = total_loss / len(test_loader)

        e_metrics(outputs_list,targets_list)
        
    return avg_loss

def e_metrics(outputs_list,targets_list):
    
    outputs_tensor = torch.cat(outputs_list, dim=0)
    targets_tensor = torch.cat(targets_list, dim=0).cpu().numpy()
    
    predicted_tensor = outputs_tensor.sigmoid() > 0.5  # 假设使用了sigmoid激活函数
    predicted_array = predicted_tensor.cpu().numpy()
    targets_array = targets_tensor
    
    metrics = {
        "accuracy": [],
        "precision": [],
        "recall": [],
        "f1": [],
        "auc": []
    }
    
    for i in range(targets_array.shape[1]):
        accuracy = accuracy_score(targets_array[:, i], predicted_array[:, i])
        precision = precision_score(targets_array[:, i], predicted_array[:, i], average=None)
        recall = recall_score(targets_array[:, i], predicted_array[:, i], average=None)
        f1 = f1_score(targets_array[:, i], predicted_array[:, i], average=None)
        auc = roc_auc_score(targets_array[:, i], outputs_tensor[:, i].cpu().numpy())
    
        metrics["accuracy"].append(accuracy)
        metrics["precision"].append(precision)
        metrics["recall"].append(recall)
        metrics["f1"].append(f1)
        metrics["auc"].append(auc)
    label_name = ['left','right']
    for i, label_metrics in enumerate(metrics["accuracy"]):
        print(f"Label {label_name[i]}:")
        for metric_name, metric_value in metrics.items():
            print(f"  {metric_name}: {metric_value[i]}")
