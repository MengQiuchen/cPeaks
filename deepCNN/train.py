# 20230412 for test all the evaluation score

import torch
import torch.nn as nn
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_auc_score

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
