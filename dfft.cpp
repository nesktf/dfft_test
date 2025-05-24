#include <fmt/core.h>

#include <complex>
#include <cstdint>
#include <span>
#include <vector>
#include <functional>

using uint32 = uint32_t;

template<typename T>
using cmplx = std::complex<T>;

template<typename T>
using span = std::span<T>;

template<typename T>
T clength(const cmplx<T>& c){
  const T sq = c.real()*c.real() + c.imag()*c.imag();
  return std::sqrt(sq);
}

template<typename T>
span<const T> to_cspan(const std::vector<T>& vec) {
  return {vec.data(), vec.size()};
}

template<typename T>
span<T> to_span(std::vector<T>& vec) {
  return {vec.data(), vec.size()};
}

template<typename F, typename T>
concept signal_fun =
  std::is_invocable_v<std::remove_cvref_t<F>, T> &&
  std::convertible_to<std::invoke_result_t<std::remove_cvref_t<F>, T>, cmplx<T>>;

template<typename T, signal_fun<T> F>
std::vector<cmplx<T>> sample_signal_n(const T& max, const T& min, uint32 n_samples, F&& signal) {
  const T dt = (max-min)/static_cast<T>(n_samples);
  std::vector<cmplx<T>> samples;
  samples.reserve(n_samples);
  for (uint32 i = 0; i < n_samples; ++i){
    const T t = min + static_cast<T>(i)*dt;
    samples.emplace_back(std::invoke(signal, t));
  }
  return samples;
}

template<typename T, typename Alloc = std::allocator<cmplx<T>>>
class dft_naive {
public:
  using value_type = T;
  using cmplx_type = cmplx<T>;

private:
  static constexpr value_type PI = static_cast<T>(M_PI);
  using cmplx_vec = std::vector<cmplx_type, Alloc>;

public:
  dft_naive(bool inverted = false)
  noexcept(std::is_nothrow_default_constructible_v<Alloc>) :
    _alloc{}, _inv{inverted} {}

  dft_naive(Alloc&& alloc, bool inverted = false)
  noexcept(std::is_nothrow_move_constructible_v<Alloc>) :
    _alloc{std::move(alloc)}, _inv{inverted} {}

  dft_naive(const Alloc& alloc, bool inverted = false)
  noexcept(std::is_nothrow_copy_constructible_v<Alloc>) :
    _alloc{alloc}, _inv{inverted} {}

public:
  template<std::random_access_iterator It>
  requires(std::same_as<typename std::iterator_traits<It>::value_type, cmplx_type>)
  void operator()(It samples_begin, It samples_end) {
    const auto n_it = std::distance(samples_begin, samples_end);
    if (n_it <= 0) {
      return;
    }

    const uint32 n = static_cast<uint32>(n_it);

    const cmplx_vec samples{samples_begin, samples_end, _alloc};
    const T omega = T{2}*PI/static_cast<T>(n) * (_inv ? T{1} : T{-1});
    const T scale = (_inv ? T{1}/static_cast<T>(n) : T{1});

    for (uint32 k = 0; k < n; ++k){
      cmplx_type sum{0, 0};
      for (uint32 i = 0; i < n; ++i) {
        const T angle = omega*static_cast<T>(k)*static_cast<T>(i);
        const cmplx_type rot{std::cos(angle), std::sin(angle)};
        sum = sum + samples[i]*rot;
      }
      *(samples_begin+k) = sum*scale;
    }
  }

public:
  const Alloc& get_allocator() const { return _alloc; }
  Alloc& get_allocator() { return _alloc; }

private:
  Alloc _alloc;
  bool _inv;
};

template<typename T, typename Alloc = std::allocator<cmplx<T>>>
class dft_ct {
public:
  using value_type = T;
  using cmplx_type = cmplx<T>;

private:
  static constexpr value_type PI = static_cast<T>(M_PI);
  using cmplx_vec = std::vector<cmplx_type, Alloc>;

public:
  dft_ct(bool inverted = false)
  noexcept(std::is_nothrow_default_constructible_v<Alloc>) :
    _alloc{}, _inv{inverted} {}

  dft_ct(Alloc&& alloc, bool inverted = false)
  noexcept(std::is_nothrow_move_constructible_v<Alloc>) :
    _alloc{std::move(alloc)}, _inv{inverted} {}

  dft_ct(const Alloc& alloc, bool inverted = false)
  noexcept(std::is_nothrow_copy_constructible_v<Alloc>) :
    _alloc{alloc}, _inv{inverted} {}

public:
  template<std::random_access_iterator It>
  requires(std::same_as<typename std::iterator_traits<It>::value_type, cmplx_type>)
  void operator()(It samples_begin, It samples_end) {
    const auto n_it = std::distance(samples_begin, samples_end);
    if (n_it <= 0) {
      return;
    }

    const uint32 n = static_cast<uint32>(n_it);
    if ((n & (n - 1))) { // Only powers of two
      return;
    }

    _compute(samples_begin, n);
  }

private:
  template<typename It>
  void _compute(It samples, uint32 n) {
    if (n == 1) {
      return;
    }

    // https://cp-algorithms.com/algebra/fft.html
    const uint32 half = n/2;
    cmplx_vec evens{_alloc};
    cmplx_vec odds{_alloc};
    evens.reserve(half);
    odds.reserve(half);
    for (uint32 i = 0; i < half; ++i) {
      evens.emplace_back(*(samples + 2*i));
      odds.emplace_back(*(samples + 2*i + 1));
    }
    _compute(odds.begin(), half);
    _compute(evens.begin(), half);

    const T angle = T{2}*PI/static_cast<T>(n) * (_inv ? T{1} : T{-1});
    const T scale = (_inv ? T{1}/T{2} : T{1});
    for (uint32 k = 0; k < half; ++k) {
      const cmplx_type w{std::cos(static_cast<T>(k)*angle), std::sin(static_cast<T>(k)*angle)};
      *(samples+k) = (evens[k] + w*odds[k])*scale;
      *(samples+k+half) = (evens[k] - w*odds[k])*scale;
    }
  }

public:
  const Alloc& get_allocator() const { return _alloc; }
  Alloc& get_allocator() { return _alloc; }

private:
  Alloc _alloc;
  bool _inv;
};

template<typename T>
class dft_inplace {
public:
  using value_type = T;
  using cmplx_type = cmplx<T>;

private:
  static constexpr value_type PI = static_cast<T>(M_PI);

public:
  dft_inplace(bool inverted = false) noexcept :
    _inv{inverted} {}

public:
  template<std::random_access_iterator It>
  requires(std::same_as<typename std::iterator_traits<It>::value_type, cmplx_type>)
  void operator()(It samples_begin, It samples_end) const {
    const auto n_it = std::distance(samples_begin, samples_end);
    if (n_it <= 0) {
      return;
    }

    const uint32 n = static_cast<uint32>(n_it);
    if ((n & (n - 1))) { // Only powers of two
      return;
    }

    // https://cp-algorithms.com/algebra/fft.html
    for (uint32 i = 1, j = 0; i < n; ++i){
      uint32 bit = n >> 1;
      for (; j & bit; bit >>= 1) {
        j ^= bit;
      }
      j ^= bit;

      if (i < j){
        std::swap(*(samples_begin+i), *(samples_begin+j));
      }
    }

    for (uint32 len = 2; len <= n; len <<= 1){
      const T angle = T{2}*PI / static_cast<T>(len) * (_inv ? T{1} : T{-1});
      const cmplx_type rot{std::cos(angle), std::sin(angle)};
      for (uint32 i = 0; i < n; i += len) {
        const uint32 half = len/2;
        cmplx_type w{T{1}};
        for (uint32 j = 0; j < half; ++j){
          const cmplx_type u = *(samples_begin+i+j);
          const cmplx_type v = (*(samples_begin+i+j+half))*w;
          *(samples_begin+i+j) = u+v;
          *(samples_begin+i+j+half) = u - v;
          w *= rot;
        }
      }
    }

    if (_inv) {
      std::for_each(samples_begin, samples_end, [&](cmplx_type& c) {
        c /= static_cast<T>(n);
      });
    }
  }

private:
  bool _inv;
};

template<typename T>
void print_samples(span<const cmplx<T>> samples, const char* msg) {
  fmt::print("{}\n", msg);
  for (uint32 i = 0; const auto& sample : samples) {
    fmt::print("- x[{}] = ({:.2f}, {:.2f}) [:.2f]\n",
               i, sample.real(), sample.imag(), clength(sample));
    ++i;
  }
  fmt::print("\n");
}

int main() {
  auto samples = sample_signal_n(2.f*M_PIf, -2.f*M_PIf, 16u, std::sin<float>);
  print_samples(to_cspan(samples), "sin(t) samples");

  dft_inplace<float> dfft;
  // dft_ct<float> dfft;
  // dft_naive<float> dfft;
  std::invoke(dfft, samples.begin(), samples.end());
  print_samples(to_cspan(samples), "sin(t) transform");
}
