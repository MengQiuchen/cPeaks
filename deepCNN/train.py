# 20230412 for test all the evaluation score

import torch
import torch.nn as nn
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm as tqdm
import numpy as np
from torch.utils.data import TensorDataset


def load_data(pos_file, neg_file,demo=True):
    # 读取正样本和负样本文件
    pos_data,neg_data = [],[]
    if demo:
        with open(pos_file, 'r') as f:
            for line in f.readlines()[:10000]:
                if not line.startswith('<'):
                    pos_data.append(line.strip())
        with open(neg_file, 'r') as f:
            for line in f.readlines()[:10000]:
                if not line.startswith('<'):
                    neg_data.append(line.strip())
    else:
        with open(pos_file, 'r') as f:
            for line in f.readlines():
                if not line.startswith('<'):
                    pos_data.append(line.strip())
        with open(neg_file, 'r') as f:
            for line in f.readlines():
                if not line.startswith('<'):
                    neg_data.append(line.strip())
        

    # 将所有序列合并为一个列表
    all_data = pos_data + neg_data

    # 初始化填充的结果列表
    padded_data = []
    label = [1]*len(pos_data)+[0]*len(neg_data)
    
    # 对每个序列进行处理
    for seq in tqdm(all_data):
        # 将碱基序列转换为数字编码，A: 0, C: 1, G: 2, T: 3, N: 4
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

        # 将onehot编码和填充后的序列添加到结果列表中
        padded_data.append(padded_seq)

    # 将结果列表转换为PyTorch张量，并返回
    
    return TensorDataset(torch.stack(padded_data),torch.tensor(np.array(label)))



def train(model, train_loader, optimizer, criterion, device):
    model.train()  # 将模型设置为训练模式
    model.to(device)
    
    total_loss = 0.0  # 总损失
    total_acc = 0.0  # 总准确率
    total_precision = 0.0  # 总精确率
    total_recall = 0.0  # 总召回率
    total_f1 = 0.0  # 总F1值
    total_auc = 0.0  # 总AUC值
    
    for batch_idx, (inputs, targets) in enumerate(train_loader):
        inputs, targets = inputs.to(device), targets.to(device)
        
        optimizer.zero_grad()  # 梯度清零
        
        outputs = model(inputs)
        loss = criterion(outputs, targets)  # 计算损失
        loss.backward()  # 反向传播
        optimizer.step()  # 更新参数
        
        _, predicted = outputs.max(1)
        acc = accuracy_score(targets.cpu().numpy(), predicted.cpu().numpy())  # 计算准确率
        precision = precision_score(targets.cpu().numpy(), predicted.cpu().numpy(), average='weighted')  # 计算精确率
        recall = recall_score(targets.cpu().numpy(), predicted.cpu().numpy(), average='weighted')  # 计算召回率
        f1 = f1_score(targets.cpu().numpy(), predicted.cpu().numpy(), average='weighted')  # 计算F1值
        auc = roc_auc_score(targets.cpu().numpy(), predicted.cpu().numpy(), average='weighted')  # 计算AUC值
        
        total_loss += loss.item()
        total_acc += acc
        total_precision += precision
        total_recall += recall
        total_f1 += f1
        total_auc += auc
        
    avg_loss = total_loss / len(train_loader)  # 平均损失
    avg_acc = total_acc / len(train_loader)  # 平均准确率
    avg_precision = total_precision / len(train_loader)  # 平均精确率
    avg_recall = total_recall / len(train_loader)  # 平均召回率
    avg_f1 = total_f1 / len(train_loader)  # 平均F1值
    avg_auc = total_auc / len(train_loader)  # 平均AUC值
    
    print('Train Loss: {:.4f}, Acc: {:.4f}, Precision: {:.4f}, Recall: {:.4f}, F1: {:.4f}, AUC: {:.4f}'.format(
        avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc))
    
    return (avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc)
    

def test(model, test_loader, criterion, device):
    model.eval()  # 将模型设置为评估模式
    model.to(device)
    
    total_loss = 0.0  # 总损失
    total_acc = 0.0  # 总准确率
    total_precision = 0.0  # 总精确率
    total_recall = 0.0  # 总召回率
    total_f1 = 0.0  # 总F1值
    total_auc = 0.0  # 总AUC值
    
    with torch.no_grad():  # 关闭梯度计算
        for batch_idx, (inputs, targets) in enumerate(test_loader):
            inputs, targets = inputs.to(device), targets.to(device)
            
            outputs = model(inputs)
            loss = criterion(outputs, targets)  # 计算损失
            
            _, predicted = outputs.max(1)
            acc = accuracy_score(targets.cpu().numpy(), predicted.cpu().numpy())  # 计算准确率
            precision = precision_score(targets.cpu().numpy(), predicted.cpu().numpy(), average='weighted')  # 计算精确率
            recall = recall_score(targets.cpu().numpy(), predicted.cpu().numpy(), average='weighted')  # 计算召回率
            f1 = f1_score(targets.cpu().numpy(), predicted.cpu().numpy(), average='weighted')  # 计算F1值
            auc = roc_auc_score(targets.cpu().numpy(), predicted.cpu().numpy(), average='weighted')  # 计算AUC值
            
            total_loss += loss.item()
            total_acc += acc
            total_precision += precision
            total_recall += recall
            total_f1 += f1
            total_auc += auc
            
    avg_loss = total_loss / len(test_loader)  # 平均损失
    avg_acc = total_acc / len(test_loader)  # 平均准确率
    avg_precision = total_precision / len(test_loader)  # 平均精确率
    avg_recall = total_recall / len(test_loader)  # 平均召回率
    avg_f1 = total_f1 / len(test_loader)  # 平均F1值
    avg_auc = total_auc / len(test_loader)  # 平均AUC值
    
    print('Test Loss: {:.4f}, Acc: {:.4f}, Precision: {:.4f}, Recall: {:.4f}, F1: {:.4f}, AUC: {:.4f}'.format(
        avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc))
    
    return (avg_loss, avg_acc, avg_precision, avg_recall, avg_f1, avg_auc)

def plot_res(train_res,test_res,num_epochs):
    
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
    plt.show()

