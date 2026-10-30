# Entegrasyon raporu: SamurayBalgami/CPPWRF → Xassemblianist/CPPWRF

Bu rapor, entegrasyonu yapacak ajan (veya geliştirici) için hazırlandı. **Bu dosya ana repoya alınmamalıdır**; yalnızca aşağıdaki 8 kod commit'i alınacak.

## Özet

| | |
|---|---|
| Kaynak | `https://github.com/SamurayBalgami/CPPWRF`, dal `ccr-165cd08e-buygym` |
| Hedef | `https://github.com/Xassemblianist/CPPWRF`, dal `main` |
| Taban (merge-base) | `fdc725c` (hedefin `main`'i, 2026-10-30 00:37 +0300) |
| Alınacak commit'ler | `fdc725c..02b872e` (8 commit, aşağıdaki tablo) |
| Alınmayacak | Bu rapor commit'i (dalın son commit'i) |
| Değişiklik | 19 dosya; kaynak kodda 2 hata düzeltmesi + 2 yorum, geri kalanı test |
| Çakışma | Yok: dal doğrudan hedefin `main`'i üzerine kurulu, fast-forward / temiz cherry-pick |

Fork'a özel `.github/workflows/sync-upstream.yml` ve ona ait 4 commit bu dala **dahil edilmedi**.

## Entegrasyon adımları

```bash
git remote add samuray https://github.com/SamurayBalgami/CPPWRF.git
git fetch samuray ccr-165cd08e-buygym
git checkout main && git pull origin main          # main hâlâ fdc725c ise
git cherry-pick fdc725c..02b872e                   # 8 commit; rapor commit'i hariç
# Alternatif (main fdc725c'de duruyorsa): git merge --ff-only 02b872e

# Doğrulama
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j && (cd build && ctest)
cmake -S . -B build64 -DCPPWRF_DOUBLE=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build64 -j && (cd build64 && ctest)
# WRF v4.8.0 -r8 derlemesi olan makinede (bu ortamda çalıştırılamadı):
cmake -S . -B build64 -DCPPWRF_DOUBLE=ON -DCPPWRF_DIFF_TESTS=ON -DWRF_DIR=/path/to/WRF && cmake --build build64 -j \
  && (cd build64 && ctest -L diff)
git push origin main
```

`main`, `fdc725c`'den ilerlemişse cherry-pick yine kullanılabilir; çakışma çıkarsa dokunulan dosyalar aşağıda listelidir.

## Commit'ler

| # | Commit | Tür | Dosyalar | Açıklama |
|---|---|---|---|---|
| 1 | `dbb8efa` | Test düzeltmesi | `tests/unit/test_core.cpp` | `test_halo_symmetric`, simetrik sınırda U/V'nin duvar yüzünü hâlâ 0 bekliyordu; `3754495`'ten beri `fill_halo` (WRF gibi) bu yüzü değiştirmiyor. **Hedefin `main`'inde `test_core` FP32 ve FP64'te başarısız** (doğrulandı: `test_core.cpp:90 xu(0, j, 0) == 0`). |
| 2 | `4b64811` | Test düzeltmesi | `tests/dyn/diffusion/driver_diffusion.cpp`, `run_test.sh` | Sürücü `defor13`'ü `ny`, `defor23`'ü `nx` adımıyla ayırıyordu; `cal_deform_and_div` ise `(.., ny)` ve `(nx, ..)` yazıyor. Sonuç: heap taşması (ASan, `deform_div.hpp:354`) ve `defor23_specified`'da 2.4e-3 fark. Artık `(nx+1)(nz+1)(ny+1)`; 130/130 alan WRF ile bit düzeyinde aynı. CUDA sürücüsü yalnızca `nvcc` ve GPU varsa çalışıyor. |
| 3 | `c23e226` | Yeni diferansiyel test | `tests/dyn/lateral_bc/*`, `CMakeLists.txt`, `src/core/grid.hpp`, `src/dyn/lateral_bc.hpp` (yalnız yorum) | `e342908` `src/bdy`'yi ve testini kaldırınca `LateralBC`'nin WRF'e karşı testi kalmamıştı. WRF'nin `relax_bdy_dry`, `spec_bdy_dry`, `spec_bdyupdate(_ph)`, `zero_grad_bdy`, `relax/spec_bdy_scalar`, `flow_dep_bdy`, `spec_bdy_final` rutinleriyle bit düzeyinde karşılaştırma; real-data ve nest; 56/56 aynı. ctest: `diff_dyn_lateral_bc`. |
| 4 | `e89ce5b` | Yeni diferansiyel test | `tests/dyn/scalars/*`, `CMakeLists.txt` | `ScalarTransport::run_stage` ↔ WRF'nin skaler taşıması (`rk_scalar_tend`, `advect_scalar(_pd)`, `rk_update_scalar(_pd)`, `set_physical_bc3d`); 4 sınır türü, tüm RK aşamaları, 5/3 ve 6/4 mertebe; 40/40 aynı, ASan/UBSan temiz. ctest: `diff_dyn_scalars`. |
| 5 | `08d2698` | Dosya izni | `tests/dyn/{lateral_bc,scalars}/run_test.sh` | `chmod +x`, diğer diferansiyel testler gibi. |
| 6 | `dc143eb` | **Hata düzeltmesi** + test | `src/dyn/diagnostics.cpp`, `tests/unit/test_diagnostics.cpp` | Bkz. Hata 2. |
| 7 | `d86e795` | Yeni diferansiyel test | `tests/dyn/halo/*`, `CMakeLists.txt` | `fill_halo`/`fill_halos` ↔ WRF'nin `set_physical_bc3d`/`bc2d`; 4 sınır türü, 3-B kütle/u/v/w ve 2-B kütle/muu/muv; 48/48 aynı. ctest: `diff_dyn_halo`. |
| 8 | `02b872e` | **Hata düzeltmesi** + test | `src/dyn/rk3.cpp`, `tests/unit/test_sound_steps.cpp` | Bkz. Hata 1. |

Her commit `(cherry picked from commit …)` satırıyla fork'taki orijinaline bağlıdır; yazar bilgisi korunmuştur.

## Kaynak koddaki hata düzeltmeleri

### Hata 1: Otomatik akustik adım sayısı WRF'nin `REAL` hassasiyetinde değil (`src/dyn/rk3.cpp`)

- WRF v4.8.0 `dyn_em/solve_em.F:462`: `num_sound_steps = max(2*(INT(300.*grid%dt/spacing - 0.01)+1), 4)`; `dt` ve `spacing` WRF `REAL`'idir (`-r8` ile double, yoksa float).
- `RK3::auto_sound_steps` bunu her zaman `double` ile hesaplıyordu. `300*dt/spacing - 0.01` tam sayıya denk geldiğinde float ve double farklı yuvarlanır:

  | dt (s) | dx (m) | float (WRF REAL(4)) | double |
  |---|---|---|---|
  | 6.7 | 1000 | 6 | 4 |
  | 20.1 | 3000 | 6 | 4 |

- Etki: FP32 derlemesi `time_step_sound = 0` ile bu durumlarda tek hassasiyetli WRF'den farklı sayıda akustik adım atıyordu. FP64 sonuçları değişmez.
- Düzeltme: hesap `real` tipinde (`adaptive_sound_steps` zaten öyleydi).

```diff
-  const double spacing = std::min(dx, dy);
-  return std::max(2 * (int(300.0 * dt / spacing - 0.01) + 1), 4);
+  const real spacing = std::min(real(dx), real(dy));
+  return std::max(2 * (int(real(300.) * real(dt) / spacing - real(0.01)) + 1), 4);
```

### Hata 2: Nemli koşumda toplam enerji `theta_m` ile hesaplanıyor (`src/dyn/diagnostics.cpp`)

- Nem varken `t_2 = theta_m - t0`. `FEnergy` sıcaklığı doğrudan bundan alıyordu, bu yüzden `c_v·T` terimi `1 + rvovrd·qv` katı kadar büyüktü (qv 20 g/kg'da ≈ %1.2).
- Düzeltme: aynı dosyadaki `FThp` gibi önce `1 + rvovrd·qv`'ye bölünüyor. `theta_mass` bilinçli olarak `theta_m` ile kalıyor (akı formundaki adveksiyonun koruduğu büyüklük).
- Yalnızca `diag.csv`'deki `total_energy` sütununu etkiler; model durumunu (prognostik alanları) etkilemez.

## Doğrulama (bu ortamda yapılanlar)

| Kontrol | Sonuç |
|---|---|
| `ctest`, FP32 (CPU, CUDA yok, netCDF yok) | 12/12 geçti |
| `ctest`, FP64 (`-DCPPWRF_DOUBLE=ON`) | 12/12 geçti |
| `test_sound_steps`, eski `rk3.cpp` ile (FP32) | 4 hata ("4 sound steps, WRF 6"), düzeltilmiş kodla geçiyor |
| `test_diagnostics`, eski `diagnostics.cpp` ile | nemli enerji hatası (1.5210e17, referans 1.5025e17), düzeltilmiş kodla geçiyor |
| `test_core`, hedefin `main`'inde | başarısız (`test_core.cpp:90`), commit 1 ile geçiyor |
| float/double yuvarlaması | Bağımsız C++ programıyla doğrulandı (yukarıdaki tablo) |
| Cherry-pick `fdc725c` üzerine | Çakışmasız; sonuç, fork `main`'inden yalnızca `sync-upstream.yml` kadar farklı |

**Doğrulanmadı:** `diff_dyn_halo`, `diff_dyn_lateral_bc`, `diff_dyn_scalars`, `diff_dyn_diffusion`. Bunlar bir WRF v4.8.0 `-r8` derlemesi (`WRF_DIR`) ve gfortran gerektiriyor, bu ortamda yok. Commit mesajlarındaki sonuçlar (56/56, 40/40, 48/48, 130/130 bit düzeyinde aynı) fork'taki önceki çalıştırmalara aittir. Birleştirmeden önce `ctest -L diff` ile tekrar çalıştırılması önerilir. WRF bulunmazsa bu testler 77 koduyla atlanır (SKIP), başarısız sayılmaz.

## Ek bulgu (bu dala dahil değil)

netCDF olmayan bir derlemede `cppwrf` programı link aşamasında hata veriyor: `undefined reference to cppwrf::cases::DomainRunner::rst_due(long) const` (`src/cases/domain_runner.cpp:710` ve devamında kullanılıyor, `:858`'de tanımlı). `src/cases` iki repoda aynı olduğu için hata hedefin `main`'inde de var. Test hedefleri etkilenmiyor. Ayrıca incelenmesi önerilir.

## Entegrasyon sonrası

Fork'un `sync-upstream.yml` iş akışı hedefin `main`'ini fork'a periyodik olarak birleştiriyor. Aynı değişiklikler iki tarafta da bulunacağı için bu birleştirmede çakışma beklenmiyor.
