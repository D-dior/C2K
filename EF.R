library(parallel)

# 'up_sigl와 dn_sigl의 모든 유전자 서명을 결합하고 중복을 제거하여 유니버설(univ) 유전자 집합 생성
up_univ = c(unlist(up_sigl) %>% unique)
dn_univ = c(unlist(dn_sigl) %>% unique)

# 두 유전자 집합(g1, g2) 및 유니버설 유전자 집합(univ)을 사용하여 EF를 계산하는 함수정의
cal_ef <- function(g1, g2, univ) {
    y = g2
    p.len = length(intersect(y, univ))
    m = intersect(g1, univ)
    D.len = length(m)
    o.len = length(intersect(y, m))
    ef = (o.len + 1) / (p.len * D.len / length(univ) + 1)
    return (ef)
}

# 병렬처리 - 1부터 10,000,000까지의 숫자 시퀀스를 생성하고, 각 숫자마다 임의의 유전자 서명을 샘플링하여 EF를 계산
ll = seq(1, 1e+7)

efl = mclapply(ll, mc.cores = 200, function(xx) {
    s1 = sample(1:length(up_sigl), 1)
    s2 = sample(1:length(up_sigl), 1)
    ef = data.frame(s1 = s1, s2 = s2,
                    up_EF = cal_ef(up_sigl[[s1]], up_sigl[[s2]], up_univ),
                    dn_EF = cal_ef(dn_sigl[[s1]], dn_sigl[[s2]], dn_univ))
    return(ef)
})
efdf = do.call(rbind, efl)
dd = efdf[, c('s1', 's2')]
efdf = efdf[!duplicated(dd), ]

# 결과 저장
saveRDS(efdf, sprintf('%s/CMAP_%s_%s_EF_res_allsample.RDS', out_dir, cell_line, dose))