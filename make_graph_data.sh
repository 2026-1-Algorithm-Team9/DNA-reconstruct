#!/bin/bash
# ============================================================
# 발표/보고서용 그래프 데이터(CSV) 자동 생성 스크립트
#
# 사용법:  bash make_graph_data.sh
#
# 하는 일:
#   - 헤더(preIndex_core.h)의 에러율/커버리지를 바꿔가며 재컴파일·실행
#   - 각 조건의 복원율(greedy, consensus)을 뽑아 CSV 2개로 저장
#   - 같은 폴더에 "에러율_복원율_(시간).csv", "커버리지_복원율_(시간).csv" 생성
#   - 원래 헤더는 끝나면 그대로 복원됨
#
# 생성된 CSV는 엑셀/구글시트에 그대로 열어 꺾은선 그래프를 만들면 됨.
# ============================================================

set -e
cd "$(dirname "$0")"

HEADER="preIndex_core.h"
BACKUP="$(mktemp)"
BIN="$(mktemp)"
cp "$HEADER" "$BACKUP"

# 각 조건을 안정적으로 보려고 REPEAT회 반복해 평균/최저를 같이 기록
REPEAT=5
# 파일명에 생성시간 포함. (한글 파일명은 일부 셸/도구에서 깨지므로 영문 사용)
STAMP=$(date "+%Y%m%d_%H%M%S")
ERR_CSV="graph_errorRate_vs_accuracy_${STAMP}.csv"
COV_CSV="graph_coverage_vs_accuracy_${STAMP}.csv"

# 헤더에서 한 줄을 값으로 치환
set_param() {  # set_param <매크로이름> <새값>
    sed -i '' "s/^#define $1 .*/#define $1     $2/" "$HEADER"
}

# 현재 헤더로 컴파일 후 1회 실행, "greedy복원율 consensus복원율"을 stdout에 출력.
# 비교요약 라인에서 콜론(:) 뒤 첫 숫자만 뽑아 공백 정렬 차이에 영향받지 않게 한다.
run_once() {
    cc -O2 -o "$BIN" preIndex.c preIndex_core.c 2>/dev/null
    out=$("$BIN" 2>/dev/null)
    g=$(echo "$out" | grep "Greedy"    | grep -oE "[0-9]+\.[0-9]+ %" | head -1 | grep -oE "[0-9]+\.[0-9]+")
    c=$(echo "$out" | grep "Consensus" | grep -oE "[0-9]+\.[0-9]+ %" | head -1 | grep -oE "[0-9]+\.[0-9]+")
    [ -z "$g" ] && g=0
    [ -z "$c" ] && c=0
    echo "$g $c"
}

# REPEAT회 실행해 greedy/consensus 각각 평균과 최저를 계산해 출력
run_repeat() {
    local gsum=0 csum=0 gmin=100 cmin=100 n=$REPEAT
    for ((i=0;i<n;i++)); do
        read g c <<< "$(run_once)"
        gsum=$(awk "BEGIN{print $gsum+$g}")
        csum=$(awk "BEGIN{print $csum+$c}")
        gmin=$(awk "BEGIN{print ($g<$gmin)?$g:$gmin}")
        cmin=$(awk "BEGIN{print ($c<$cmin)?$c:$cmin}")
        sleep 1.05   # srand(time)이 1초 단위라 매번 다른 데이터를 위해
    done
    gavg=$(awk "BEGIN{printf \"%.2f\", $gsum/$n}")
    cavg=$(awk "BEGIN{printf \"%.2f\", $csum/$n}")
    echo "$gavg $gmin $cavg $cmin"
}

echo "데이터 수집 시작 (조건당 ${REPEAT}회 반복)..."

# ---------- 그래프 A: 에러율 vs 복원율 (커버리지는 기본값 고정) ----------
echo "[1/2] 에러율 실험..."
echo "error_rate_percent,greedy_avg,greedy_min,consensus_avg,consensus_min" > "$ERR_CSV"
for E in 0 1 2 3 5 8; do
    set_param "ERROR_RATE_PERCENT" "$E"
    read gavg gmin cavg cmin <<< "$(run_repeat)"
    echo "$E,$gavg,$gmin,$cavg,$cmin" >> "$ERR_CSV"
    echo "  에러율 ${E}% -> greedy ${gavg}%, consensus ${cavg}%"
done
cp "$BACKUP" "$HEADER"   # 에러율 원복

# ---------- 그래프 B: 커버리지 vs 복원율 (에러율 1% 고정) ----------
echo "[2/2] 커버리지 실험 (에러율 1% 고정)..."
set_param "ERROR_RATE_PERCENT" "1"
echo "coverage,greedy_avg,greedy_min,consensus_avg,consensus_min" > "$COV_CSV"
for COV in 20 30 40 60 80 100; do
    set_param "COVERAGE" "$COV"
    read gavg gmin cavg cmin <<< "$(run_repeat)"
    echo "$COV,$gavg,$gmin,$cavg,$cmin" >> "$COV_CSV"
    echo "  커버리지 ${COV}배 -> consensus 평균 ${cavg}%, 최저 ${cmin}%"
done

# 헤더 원복 + 임시파일 정리
cp "$BACKUP" "$HEADER"
rm -f "$BACKUP" "$BIN"

echo ""
echo "완료! 생성된 파일:"
echo "  - $ERR_CSV   (에러율 vs 복원율 그래프용)"
echo "  - $COV_CSV   (커버리지 vs 복원율 그래프용)"
echo "엑셀/구글시트로 열어 꺾은선 그래프를 만드세요."
