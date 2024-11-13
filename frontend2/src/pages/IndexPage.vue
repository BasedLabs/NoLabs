<template>
  <q-page class="row items-center justify-evenly">
    <example-component
      title="Example component"
      active
      :todos="todos"
      :meta="meta"
    ></example-component>
  </q-page>
</template>

<script lang="ts">
import { defineComponent } from 'vue';
import { Todo, Meta } from 'components/models';
import ExampleComponent from 'components/ExampleComponent.vue';

export default defineComponent({
  name: 'IndexPage',

  components: {
    ExampleComponent
  },

  data () {
    const todos: Todo[] = [
      {
        id: 1,
        content: 'ct1'
      },
      {
        id: 2,
        content: 'ct2'
      },
      {
        id: 3,
        content: 'ct3'
      },
      {
        id: 4,
        content: 'ct4'
      },
      {
        id: 5,
        content: 'ct5'
      }
    ];

    const meta: Meta = {
      totalCount: 1200
    };

    return { todos, meta };
  }
});
</script>
