#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define refLength 100   // 원본 염기서열 길이
#define fragLength 10   // 조각 길이
#define fragNum 10      // 조각 개수

char* makeRef() {
    char basis[4] = { 'A', 'T', 'C', 'G' };

    // 염기 서열 저장하기 위한 배열 동적 할당 (+1은 문자열 종료 문자 '\0'를 위해)
    char* ref = (char*)malloc((refLength + 1) * sizeof(char));
    srand(time(NULL));

    for (int i = 0; i < refLength; i++) {
        int index = rand() % 4;
        ref[i] = basis[index];
    }
    ref[refLength] = '\0'; // 문자열 종료 문자 추가

    return ref;
}

char** makeFrag() {
    // A, T, C, G 염기를 배열로 정의
    char bases[4] = {'A', 'T', 'C', 'G'};
    
    // 조각들을 저장하기 위한 2차원 배열 동적 할당
    char** frags = (char**)malloc(fragNum * sizeof(char*)); 

    for (int i = 0; i < fragNum; i++) {
        // 각 조각을 저장하기 위한 배열 동적 할당 (+1은 '\0' 자리를 위함)
        frags[i] = (char*)malloc((fragLength + 1) * sizeof(char)); 

        for (int j = 0; j < fragLength; j++) {
            // 0, 1, 2, 3 중 랜덤한 인덱스를 뽑아 A, T, C, G 중 하나를 무작위 선택
            int index = rand() % 4; 
            frags[i][j] = bases[index]; 
        }
        frags[i][fragLength] = '\0'; // 문자열 종료 문자 추가
    }

    return frags;
}

// 문자를 0~3 숫자로 변환 (A, C, G, T 기반 2비트 인코딩)
int char_to_int(char c) {
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'T') return 3;
    return 0; // 예외 처리용 기본값
}

// k-mer에 대한 연속 해시값 계산
int get_hash(const char* str, int k) {
    int hash = 0;
    for (int i = 0; i < k; i++) {
        hash = (hash << 2) | char_to_int(str[i]);
    }
    return hash;
}

// 기본 해싱 알고리즘 (Rolling Hash + Prefix Sum 기반 정적 배열)
void hashing_search(const char* genome, const char* read, int k, int threshold) {
    int g_len = strlen(genome);
    int r_len = strlen(read);

    // k-mer 길이가 read 길이보다 길면 탐색 불가
    if (k > r_len) return;

    // 해시 크기 계산 (4^k)
    int hash_size = 1 << (2 * k);

    // 빈도수 저장을 위한 배열 할당 및 초기화 (0으로 세팅)
    int* count = (int*)calloc(hash_size + 1, sizeof(int));

    // 단계 1: k-mer 빈도수 카운팅
    for (int i = 0; i <= g_len - k; i++) {
        int h = get_hash(genome + i, k);
        count[h]++;
    }

    // 단계 2: 누적합 (Prefix Sum) 계산을 통한 원본 인덱스 저장 좌표 산출
    int* prefix_sum = (int*)calloc(hash_size + 1, sizeof(int));

    for (int i = 1; i <= hash_size; i++) {
        prefix_sum[i] = prefix_sum[i - 1] + count[i - 1];
    }

    // 단계 3: 정적 1차원 배열 할당 및 실제 위치(인덱스) 기록
    int total_kmers = g_len - k + 1;
    int* index_array = (int*)malloc(total_kmers * sizeof(int));
    int* current_offset = (int*)calloc(hash_size, sizeof(int)); // 위치 중복 방지용 오프셋

    for (int i = 0; i <= g_len - k; i++) {
        int h = get_hash(genome + i, k);
        int pos = prefix_sum[h] + current_offset[h];
        index_array[pos] = i;
        current_offset[h]++;
    }

    // 단계 4: 탐색 수행 (read의 첫 번째 k-mer를 시드로 사용)
    int read_hash = get_hash(read, k);
    int start = prefix_sum[read_hash]; // 해당 해시값이 모여있는 시작점
    int count_match = count[read_hash]; // 해당 해시값의 총 개수

    // 일치하는 시드 위치에서만 미스매치 검사 진행
    for (int i = 0; i < count_match; i++) {
        int g_idx = index_array[start + i];

        // 탐색 범위가 원본 길이를 초과하면 무시
        if (g_idx + r_len > g_len) continue;

        int mismatch = 0;

        for (int j = 0; j < r_len; j++) {
            if (genome[g_idx + j] != read[j]) {
                mismatch++;
                if (mismatch > threshold) break;
            }
        }

        // 최종 조건 통과 시 출력
        if (mismatch <= threshold) {
            printf("Match at index %d (mismatches: %d)\n", g_idx, mismatch);
        }
    }

    // 동적 할당 메모리 해제 (누수 방지)
    free(count);
    free(prefix_sum);
    free(index_array);
    free(current_offset);
}

int main() {

    // 테스트용 시퀀스 무작위 생성
    char* genome = makeRef();
    char** reads = makeFrag(genome);

    int mismatch_threshold = 1; // 1개까지 불일치 허용
    int k_mer_size = 4;         // 부분 서열(k-mer) 길이는 4로 설정

    printf("원본 서열: %s\n\n", genome);

    // 10개의 조각에 대해 탐색 실행
    for (int i = 0; i < fragNum; i++) {
        printf("--- [%d번째 조각: %s] Hashing Search ---\n", i + 1, reads[i]);
        hashing_search(genome, reads[i], k_mer_size, mismatch_threshold);
        printf("\n");
    }

    // 메모리 해제
    for (int i = 0; i < fragNum; i++) {
        free(reads[i]); // 개별 조각 해제
    }
    free(reads);  // 조각 포인터 배열 해제
    free(genome); // 원본 서열 해제

    return 0;

}