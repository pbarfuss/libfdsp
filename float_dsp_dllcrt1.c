#include <windef.h>
#include <winbase.h>
#include <winnt.h>
LONG NTAPI NtQueryPerformanceCounter(LARGE_INTEGER *PerformanceCounter, LARGE_INTEGER *PerformanceFrequency);
LONG NTAPI LdrDisableThreadCalloutsForDll(VOID *BaseAddress);

#ifdef _WIN64
#define DEFAULT_SECURITY_COOKIE 0x00002B992DDFA232ll
#else
#define DEFAULT_SECURITY_COOKIE 0xBB40E64E
#endif
DECLSPEC_SELECTANY UINT_PTR __security_cookie = DEFAULT_SECURITY_COOKIE;
DECLSPEC_SELECTANY UINT_PTR __security_cookie_complement = ~(DEFAULT_SECURITY_COOKIE);

static void __security_init_cookie (void)
{
  UINT_PTR cookie;
  LARGE_INTEGER perfctr, perffreq;

  if (__security_cookie != DEFAULT_SECURITY_COOKIE) {
      __security_cookie_complement = ~__security_cookie;
      return;
  }

  NtQueryPerformanceCounter(&perfctr, &perffreq);                                                                                                                                                                 
#ifdef _WIN64                                                                                                                                                                                                     
  cookie = perfctr.QuadPart;
#else
  cookie = perfctr.LowPart;
  cookie ^= perfctr.HighPart;
#endif

  if (cookie == DEFAULT_SECURITY_COOKIE)
    cookie = DEFAULT_SECURITY_COOKIE + 1;
  __security_cookie = cookie;
  __security_cookie_complement = ~cookie;
}

WINBOOL WINAPI
DllMainCRTStartup (HANDLE hDllHandle, DWORD dwReason, LPVOID lpreserved)
{
  if (dwReason == DLL_PROCESS_ATTACH) {
    __security_init_cookie();
    LdrDisableThreadCalloutsForDll(hDllHandle);
  }
  return TRUE;
}

