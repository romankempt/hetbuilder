<script lang="ts">
	export let name: string = "Roman";

	let rand = -1;
	function getRand() {
		fetch("./rand")
			.then((d) => d.text())
			.then((d) => (rand = <number>d));
	}

	let postVar;
	let fileVar;

	function submitForm() {
		event.preventDefault();

		const dataArray = new FormData();
		dataArray.append("superHeroName", postVar);
		dataArray.append("uploadFile", fileVar);

		fetch("./upload", {
			method: "POST",
			headers: [["Content-Type", "multipart/form-data"]],
			body: dataArray,
		})
			.then((response) => {
				console.log(response);
			})
			.catch((error) => {
				// Upload failed
			});
	}
</script>

<main>
	<div>
		<h1>Hello {name} from me!</h1>
		<p>
			Visit the <a href="https://svelte.dev/tutorial">Svelte tutorial</a> to
			learn how to build Svelte apps.
		</p>
	</div>

	<div>
		<p>
			Your number is this one here {rand}!
			<button on:click={getRand}>Get a random number</button>
		</p>
	</div>

	<div>
		<form on:submit={submitForm}>
			<input
				type="text"
				bind:value={postVar}
				placeholder={"Superhero Name"}
			/>
			<br />
			<input type="file" bind:files={fileVar} />
			<br />
			<input type="submit" />
		</form>
	</div>
</main>

<style>
	main {
		text-align: center;
		padding: 1em;
		max-width: 240px;
		margin: 0 auto;
	}

	h1 {
		color: #ff3e00;
		text-transform: uppercase;
		font-size: 4em;
		font-weight: 100;
	}

	button {
		color: "blue";
	}

	@media (min-width: 640px) {
		main {
			max-width: none;
		}
	}
</style>
