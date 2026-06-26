# Mini Time-Domain Analysis（迷你时域分析）

**从零开始、零依赖的 C 语言实现**，涵盖经典控制理论的时域分析。每个模块对应 MIT（及其他顶尖大学）的一门或多门课程，将教科书中的公式转化为可运行的 C 代码，实现理论与实践的桥接。

## 子模块总览

| 子模块 | 主题 | 参考课程 |
|--------|--------|-------------|
| [mini-higher-order-reduction](mini-higher-order-reduction/) | 平衡截断、主导极点分析、模态截断、模型降阶 | MIT 6.242, ETH 227-0216 |
| [mini-impulse-step-response](mini-impulse-step-response/) | 卷积积分、离散时间响应、LTI 系统冲激/阶跃响应 | MIT 6.241J, Nise 第4章 |
| [mini-second-order-systems](mini-second-order-systems/) | 标准二阶系统、ζ/ωn 参数映射、解析与数值响应 | MIT 6.241J, Nise 第4-5章 |
| [mini-sensitivity-analysis](mini-sensitivity-analysis/) | Bode 灵敏度积分、特征值/参数/时域/传递函数灵敏度 | MIT 6.241J, Skogestad, Caltech CDS 110 |
| [mini-stability-routh-hurwitz](mini-stability-routh-hurwitz/) | Routh-Hurwitz 判据、Jury 稳定性（离散时间）、相对稳定裕度 | MIT 6.241J, Nise 第6章, Stanford ENGR 105 |
| [mini-steady-state-error](mini-steady-state-error/) | 误差常数（Kp/Kv/Ka）、系统型别、扰动抑制、终值定理 | MIT 6.241J, Nise 第7章, Stanford ENGR 105 |
| [mini-time-domain-design](mini-time-domain-design/) | 特征值配置、状态空间设计、矩阵函数、线性代数运算 | MIT 6.241J, Nise 第9章, Stanford ENGR 105 |
| [mini-transient-specs](mini-transient-specs/) | 上升/调节时间、超调量、一阶/二阶/高阶瞬态、ISE/IAE/ITAE 指标 | MIT 6.241J, Nise 第4章, Stanford ENGR 105 |

## 设计理念

- **零外部依赖** — 纯 C（C99/C11），仅使用 `libc` 和 `libm`
- **模块自包含** — 每个目录自带 `include/`、`src/`、`examples/`、`tests/`
- **理论到代码的映射** — 每个模块包含 `docs/` 目录，内有课程对齐说明及教科书参考文献（Nise、Ogata、Franklin）
- **经典控制传承** — 覆盖控制教育的时域支柱：瞬态分析、误差分析、稳定性与灵敏度

## 构建方式

每个模块相互独立。进入模块目录后运行：

```bash
cd mini-transient-specs
make all    # 构建全部
make test   # 运行测试
```

需要 **GCC** 和 **GNU Make**。

## 项目结构

```
mini-time-domain-analysis/
├── mini-higher-order-reduction/   # 平衡截断、主导极点、模态模型降阶
├── mini-impulse-step-response/    # 卷积积分、冲激/阶跃响应、离散时间响应
├── mini-second-order-systems/     # 标准二阶 LTI 系统、ζ/ωn 参数化
├── mini-sensitivity-analysis/     # Bode 积分、特征值/参数/时域灵敏度分析
├── mini-stability-routh-hurwitz/  # Routh-Hurwitz、Jury 判据、相对稳定裕度
├── mini-steady-state-error/       # 误差常数、系统型别、扰动抑制、终值定理
├── mini-time-domain-design/       # 特征值设计、状态空间方法、矩阵函数
└── mini-transient-specs/          # 瞬态响应指标、性能指数（ISE/IAE/ITAE）
```

## 许可证

MIT
