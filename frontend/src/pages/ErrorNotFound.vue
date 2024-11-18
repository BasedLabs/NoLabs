<!-- src/views/NotFound.vue -->
<template>
  <div class="not-found">
    <h1>404</h1>
    <p>Oops! The page you’re looking for doesn’t exist.</p>
    <router-link to="/">Go back home</router-link>
  </div>
</template>

<script>
export default {
  name: 'NotFound',
};
</script>

<style scoped>
.not-found {
  text-align: center;
  padding: 50px;
}

.not-found h1 {
  font-size: 96px;
  margin-bottom: 0;
}

.not-found p {
  font-size: 24px;
  margin: 10px 0 20px;
}

.not-found a {
  color: #42b983;
  font-size: 18px;
  text-decoration: underline;
}
</style>
