#!/bin/bash
# ============================================================
# 탐색(매핑) 벤치마크 데이터(CSV) 자동 생성 스크립트
#   대상: Trivial 완전탐색  vs  Simple Hash 인덱스 매핑
#
# 사용법:  bash make_search_benchmark.sh
#
# 하는 일:
#   - 두 파일의 에러율(ERROR_RATE_PERCENT)을 0~5%로 바꿔가며 재컴파일·실행
#   - 각 조건의 (탐색 시간 / 위치 정확도)를 뽑아 CSV 2개로 저장
#   - 같은 폴더에 search_speed_(시간).csv, search_accuracy_(시간).csv 생성
#   - 원래 소스(.c)는 끝나면 그대로 복원됨
#
# 생성 CSV는 엑셀/구글시트로 열어 그래프(속도 비교 / 정확도 비교)를 만들면 됨.
# ============================================================

set -e
cd "$(dirname "$0")"

TRI="trivial_search.c"
HSH="simple_hashing_search.c"
TRI_BAK="$(mktemp)"; HSH_BAK="$(mktemp)"
BIN="$(mktemp)"
cp "$TRI" "$TRI_BAK"; cp "$HSH" "$HSH_BAK"

# 어떤 이유로 중단돼도 원본 소스를 반드시 복원 (안전장치)
restore() { cp "$TRI_BAK" "$TRI"; cp "$HSH_BAK" "$HSH"; rm -f "$TRI_BAK" "$HSH_BAK" "$BIN"; }
trap restore EXIT

REPEAT=3                              # 조건당 반복 횟수(평균). trivial이 느려 과하지 않게.
STAMP=$(date "+%Y%m%d_%H%M%S")
SPEED_CSV="search_speed_${STAMP}.csv"
ACC_CSV="search_accuracy_${STAMP}.csv"

# 한 소스의 ERROR_RATE_PERCENT 값을 바꾼다
set_err() {  # set_err <file> <value>
    sed -i '' "s/^#define ERROR_RATE_PERCENT [0-9]*/#define ERROR_RATE_PERCENT $2/" "$1"
}

# 컴파일 후 1회 실행 → "시간 위치정확도" 출력
run_one() {  # run_one <srcfile>
    cc -O2 -o "$BIN" "$1" 2>/dev/null
    out=$("$BIN" 2>/dev/null)
    t=$(echo "$out" | grep "탐색 시간" | grep -oE "[0-9]+\.[0-9]+" | head -1)
    a=$(echo "$out" | grep "위치 정확도" | grep -oE "[0-9]+\.[0-9]+" | head -1)
    [ -z "$t" ] && t=0; [ -z "$a" ] && a=0
    echo "$t $a"
}

# REPEAT회 평균 → "시간평균 정확도평균"  (awk -v로 안전하게 계산)
run_avg() {  # run_avg <srcfile>
    local ts=0 as=0 n=$REPEAT t a
    for ((i=0;i<n;i++)); do
        read t a <<< "$(run_one "$1")"
        ts=$(awk -v x="$ts" -v y="$t" 'BEGIN{printf "%.6f", x+y}')
        as=$(awk -v x="$as" -v y="$a" 'BEGIN{printf "%.4f", x+y}')
        sleep 1.05      # srand(time)이 1초 단위 → 매번 다른 데이터
    done
    local tavg aavg
    tavg=$(awk -v s="$ts" -v n="$n" 'BEGIN{printf "%.6f", s/n}')
    aavg=$(awk -v s="$as" -v n="$n" 'BEGIN{printf "%.2f", s/n}')
    echo "$tavg $aavg"
}

echo "탐색 벤치마크 데이터 수집 시작 (조건당 ${REPEAT}회 반복)..."
echo "error_rate,trivial_time,simple_time,speedup" > "$SPEED_CSV"
echo "error_rate,trivial_accuracy,simple_accuracy" > "$ACC_CSV"

for E in 0 1 2 3 5; do
    set_err "$TRI" "$E"; set_err "$HSH" "$E"
    read tt ta <<< "$(run_avg "$TRI")"      # trivial: 시간, 정확도
    read st sa <<< "$(run_avg "$HSH")"      # simple : 시간, 정확도
    speedup=$(awk -v a="$tt" -v b="$st" 'BEGIN{ if (b>0) printf "%.1f", a/b; else print "NA" }')
    echo "$E,$tt,$st,$speedup" >> "$SPEED_CSV"
    echo "$E,$ta,$sa" >> "$ACC_CSV"
    echo "  에러 ${E}% -> trivial ${tt}s / simple ${st}s (${speedup}x), 정확도 tri ${ta}% / hash ${sa}%"
done

# 소스 원복은 trap(restore)이 종료 시 자동 수행

echo ""
echo "완료! 생성된 파일:"
echo "  - $SPEED_CSV     (에러율별 속도: trivial vs simple + 배속)"
echo "  - $ACC_CSV       (에러율별 위치 정확도: trivial vs simple)"
echo "엑셀/구글시트로 열어 그래프를 만드세요."
