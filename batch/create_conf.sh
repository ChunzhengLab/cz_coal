#!/bin/bash
# 初始化临时文件来记录每个日期文件夹的.conf文件数量
conf_counts_file=$(mktemp)
seq_numbers_file=$(mktemp)
# 清空或创建submit.sh文件
: >submit.sh
# 初始化.conf文件计数器
conf_count=0
# 读取 zpc.list 文件
while IFS= read -r line; do
  # 使用正则表达式提取日期和序号信息
  if [[ $line =~ ([0-9]{10})/zpc-([0-9]+)\.root ]]; then
    date_folder="${BASH_REMATCH[1]}"
    seq_number="${BASH_REMATCH[2]}"
    # 确保记录中有此日期文件夹的条目
    if ! grep -q "^$date_folder " "$conf_counts_file"; then
      echo "$date_folder 0" >>"$conf_counts_file"
      echo "$date_folder " >>"$seq_numbers_file"
    fi
    # 创建并进入目标输出目录
    output_dir="/storage/fdunphome/wangchunzheng/cz_coal/batch/${date_folder}"
    mkdir -p "$output_dir"
    cd "$output_dir" || exit 1
    # 创建 .conf 文件
    conf_file="${seq_number}.conf"
    cat >"$conf_file" <<EOL
isDebug=false
isWriteEvents=true
isCalculateObvs=true
isRemoveHFQuarks=true
eventType=kAMPT
r_bm=0.5
flavourBreakTolerance=0.0
coalescenceAlgorithm=kFromParton
inputFile=${line}
outputFile=${output_dir}/data_coalHadrons-${seq_number}.root
obvsFile=${output_dir}/obvs_coalHadrons-${seq_number}.root
EOL
    echo "Created config file: $conf_file in $output_dir"
    # 增加计数器
    awk -v df="$date_folder" '{if ($1 == df) $2 += 1}1' "$conf_counts_file" >tmpfile && mv tmpfile "$conf_counts_file"
    # 记录文件序号
    awk -v df="$date_folder" -v seq="$seq_number" '{if ($1 == df) $2 = $2 " " seq}1' "$seq_numbers_file" >tmpfile && mv tmpfile "$seq_numbers_file"
    # 增加 .conf 文件计数
    ((conf_count++))
    # 返回上一级目录
    cd - >/dev/null || exit 1
  else
    echo "No match for line: $line"
  fi
done <zpc_sorted.list
# 创建 condor.sh 和 condor.con 文件
while read -r date_folder count; do
  # 获取该日期文件夹的所有序号并排序
  seq_numbers=$(awk -v df="$date_folder" '$1 == df {print $2}' "$seq_numbers_file")
  IFS=' ' read -r -a seq_array <<<"$seq_numbers"
  sorted_seq=$(echo "${seq_array[@]}" | tr ' ' '\n' | sort -n)
  # 获取最大序号
  max_seq=$(echo "$sorted_seq" | tail -n 1)
  # 创建 condor.sh 文件
  output_dir="/storage/fdunphome/wangchunzheng/cz_coal/batch/${date_folder}"
  cat >"${output_dir}/condor.sh" <<EOL
#!/bin/bash
# Source the ROOT environment
source /opt/root61404/bin/thisroot.sh
# Ensure an index is provided as an argument
if [ -z "\$1" ]; then
  echo "Usage: \$0 <index>"
  exit 1
fi
Index=\$1
# Config file should be in the current directory and named as <index>.conf
Configuration="./\${Index}.conf"
# Path to the executable
Execute="../../bin/Coalescence"
# Check if the configuration file exists
if [ ! -f "\$Configuration" ]; then
  echo "Configuration file \${Configuration} not found!"
  exit 1
fi
# Check if the executable exists and is executable
if [ ! -x "\$Execute" ]; then
  echo "Executable \${Execute} not found or not executable!"
  exit 1
fi
\$Execute "\$Configuration"
EOL
  # 设置 condor.sh 文件权限为 755
  chmod 755 "${output_dir}/condor.sh"
  echo "Created and set permissions for condor.sh in $output_dir"
  # 创建 condor.con 文件
  cat >"${output_dir}/condor.con" <<EOL
Universe     = vanilla
Notification = Error
Initialdir   = ${output_dir}
Executable   = \$(Initialdir)/condor.sh
Arguments    = \$(Process)
Log          = \$(Initialdir)/\$(Process).log
Output       = \$(Initialdir)/\$(Process).out
Error        = \$(Initialdir)/\$(Process).err
Notify_user  = wangcz22@m.fudan.edu.cn
Queue $max_seq
EOL
  echo "Created condor.con in $output_dir"
  # 添加提交命令到submit.sh
  echo "condor_submit ${output_dir}/condor.con" >>submit.sh
done <"$conf_counts_file"
# 设置 submit.sh 文件权限为 755
chmod 755 submit.sh
# 输出统计信息
echo "Total number of .conf files created: $conf_count"
# 删除临时文件
rm "$conf_counts_file" "$seq_numbers_file"