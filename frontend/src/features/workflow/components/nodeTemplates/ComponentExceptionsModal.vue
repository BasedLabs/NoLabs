<script lang="ts">
import {defineComponent} from 'vue'
import {QSpinnerOrbit} from "quasar";

export default defineComponent({
  name: "ComponentExceptionsModal",
  props: {
    exceptions: {
      type: Array<string>,
      required: true
    },
    visible: {
      type: Boolean,
      required: true
    }
  },
  emits: ['update:visible'],
  computed: {
    show: {
      get() {
        return this.visible;
      },
      set(val: boolean) {
        this.$emit('update:visible', val)
      }
    },
    exceptionsString(): string {
      if(!this.exceptions || this.exceptions.length === 0) {
        return '';
      }

      return this.exceptions.join('\n-------------------------------------------------\n');
    },
  },
  methods: {
    onClose() {
      this.$emit('update:visible', false);
    }
  }
})
</script>

<template>
  <q-dialog
    v-model="show"
    auto-close
  >
    <q-card style="width: 90vh; max-width: 95vw; height: 90vh; max-height: 90vh">
      <q-card-section>
        <div class="text-h6">Last exceptions</div>
      </q-card-section>

      <q-card-section class="q-pt-none">
        <div class="row q-col-gutter-xs" v-if="exceptions">
          <div class="col-xs-12 col-md-12">
            <q-input
              v-model="exceptionsString"
              readonly
              filled
              autogrow
            />
          </div>
        </div>
      </q-card-section>
    </q-card>
  </q-dialog>

</template>


