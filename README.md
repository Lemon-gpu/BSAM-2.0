> TLDR: [BSAM](https://deepwiki.com/ZhenlinGuo/BSAM-2.0)
> 
> Author of this note: ChatGPT 5, Gemini Pro 2.5, VSCode Copilot, Huan Chen

# 阅读代码的笔记

## 10-16-2025：
    Core: `bsamroutines.f90` 中的`AMR()`直接体现了自适应密化的核心逻辑，因为有递归
    Todo: 接着看AMR()

## 10-17-2025：
    Core: 数据结构不是一个标准的n叉树，它包含与V-Cycle相关的“负层”。网格被序列化处理，每一层都组织成一个链表。
    Todo: 继续研究其具体的数据结构实现。

## 10-18-2025 & 10-19-2025：
    Core: 🐟
    
