<script lang="ts">
import { defineComponent } from 'vue';
import ExperimentsView from './ExperimentsView.vue';

export default defineComponent({
    name: 'ProteinDesignPage',
    mounted() {
    },
    components: {
        ExperimentsView
    }
});
</script>

<template>
    <ExperimentsView />
</template>