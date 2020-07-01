#ifdef TRACE_EXTRAE
#include "extrae_user_events.h"
#include "extrae_types.h"

#define EVENT_NONE 0
#define EVENT_ALLOCATE  1
#define EVENT_INIT 2
#define EVENT_CREATE 3
#define EVENT_COMPUTE 4
#define EVENT_SIMULATION 5
#define EVENT_VALIDATION 6

static extrae_type_t et;
static extrae_value_t ev[7] = {0, 10, 11, 12, 13, 14, 15};
static char *extrae_names[] = {"none", "matrix_allocate", "matrix_init", "matrix_create_tasks", "matrix_compute", "complete_simulation", "matrix_validate"};

#define REGISTER_EXTRAE() do { \
  unsigned nvalues = 7; \
  Extrae_define_event_type(&et, "Operations", &nvalues, ev, extrae_names); \
} while(0)
#define EXTRAE_ENTER(_e) Extrae_eventandcounters(et, ev[_e])
#define EXTRAE_EXIT(_e)  Extrae_eventandcounters(et, ev[EVENT_NONE])
#else
#define EXTRAE_ENTER(_e)
#define EXTRAE_EXIT(_e)
#define REGISTER_EXTRAE()
#endif