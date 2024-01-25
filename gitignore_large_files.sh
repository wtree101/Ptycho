#!/bin/bash

# 定义大文件的大小阈值（以字节为单位）
large_file_threshold=1000000

# 检测项目文件夹中的大文件
echo "检测大文件..."
large_files=$(find . -type f -size +${large_file_threshold}c)

if [ -z "$large_files" ]; then
  echo "没有找到大文件。"
else
  echo "找到以下大文件："
  echo "$large_files"
  
  # 将大文件添加到.gitignore中
  echo "将大文件添加到.gitignore中..."
  echo "$large_files" >> .gitignore
  
  # 从Git缓存中删除大文件
  echo "从Git缓存中删除大文件..."
  git rm --cached $large_files
  
  echo "完成！大文件已被忽略并从Git缓存中删除。"
fi